def read3d(path):
    arr = tiff.imread(str(path))
    arr = np.asarray(arr)

    if arr.ndim == 2:
        arr = arr[np.newaxis, :, :]

    if arr.ndim == 4 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]

    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D volume, got shape={arr.shape} for {path}")

    return arr


def label_counts(vol):
    ids, counts = np.unique(vol, return_counts=True)
    mask = ids > 0
    return ids[mask].astype(np.int64), counts[mask].astype(np.int64)


def border_label_set(vol):
    vals = np.concatenate([
        vol[0, :, :].ravel(),
        vol[-1, :, :].ravel(),
        vol[:, 0, :].ravel(),
        vol[:, -1, :].ravel(),
        vol[:, :, 0].ravel(),
        vol[:, :, -1].ravel(),
    ])
    vals = vals[vals > 0]
    return set(vals.astype(np.int64).tolist())


def write_csv(path, header, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gt", required=True, help="Ground-truth label crop")
    ap.add_argument("--pred", required=True, help="Predicted label crop")
    ap.add_argument("--markers", default=None, help="Optional marker components")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-gt-area", type=int, default=1)
    ap.add_argument("--ignore-gt-border", action="store_true")
    ap.add_argument("--iou-thresholds", default="0.10,0.25,0.50")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gt = read3d(args.gt).astype(np.int64)
    pred = read3d(args.pred).astype(np.int64)

    if gt.shape != pred.shape:
        raise RuntimeError(f"Shape mismatch: gt={gt.shape}, pred={pred.shape}")

    print(f"[eval] gt shape   = {gt.shape}, dtype={gt.dtype}")
    print(f"[eval] pred shape = {pred.shape}, dtype={pred.dtype}")

    gt_ids_all, gt_counts_all = label_counts(gt)
    pred_ids, pred_counts = label_counts(pred)

    border_ids = border_label_set(gt)

    valid_gt = []
    gt_count_dict = {}
    for gid, cnt in zip(gt_ids_all, gt_counts_all):
        if cnt < args.min_gt_area:
            continue
        if args.ignore_gt_border and int(gid) in border_ids:
            continue
        valid_gt.append(int(gid))
        gt_count_dict[int(gid)] = int(cnt)

    valid_gt = np.array(valid_gt, dtype=np.int64)

    pred_count_dict = {int(pid): int(cnt) for pid, cnt in zip(pred_ids, pred_counts)}

    gt_fg = gt > 0
    pred_fg = pred > 0

    inter_fg = int(np.logical_and(gt_fg, pred_fg).sum())
    union_fg = int(np.logical_or(gt_fg, pred_fg).sum())
    gt_fg_count = int(gt_fg.sum())
    pred_fg_count = int(pred_fg.sum())

    fg_iou = inter_fg / union_fg if union_fg else 0.0
    fg_dice = (2 * inter_fg) / (gt_fg_count + pred_fg_count) if (gt_fg_count + pred_fg_count) else 0.0
    fg_precision = inter_fg / pred_fg_count if pred_fg_count else 0.0
    fg_recall = inter_fg / gt_fg_count if gt_fg_count else 0.0

    # Pair intersections between valid GT labels and predicted labels.
    valid_gt_set = set(valid_gt.tolist())
    pair_mask = (gt > 0) & (pred > 0)

    if valid_gt.size > 0:
        pair_mask &= np.isin(gt, valid_gt)

    gt_pair = gt[pair_mask]
    pred_pair = pred[pair_mask]

    pair_rows = []
    intersections = {}

    if gt_pair.size > 0:
        max_pred = int(pred.max()) + 1
        codes = gt_pair.astype(np.int64) * max_pred + pred_pair.astype(np.int64)
        unique_codes, inter_counts = np.unique(codes, return_counts=True)

        for code, inter in zip(unique_codes, inter_counts):
            gid = int(code // max_pred)
            pid = int(code % max_pred)
            inter = int(inter)
            intersections[(gid, pid)] = inter
            pair_rows.append((gid, pid, inter))

    write_csv(
        out_dir / "pair_intersections.csv",
        ["gt_id", "pred_id", "intersection_voxels"],
        pair_rows,
    )

    # Best prediction for each GT.
    gt_best_rows = []
    matched_gt_at = {0.10: 0, 0.25: 0, 0.50: 0}

    by_gt = {}
    by_pred = {}

    for (gid, pid), inter in intersections.items():
        by_gt.setdefault(gid, []).append((pid, inter))
        by_pred.setdefault(pid, []).append((gid, inter))

    for gid in valid_gt:
        gid = int(gid)
        gcount = gt_count_dict.get(gid, int((gt == gid).sum()))
        candidates = by_gt.get(gid, [])

        if not candidates:
            gt_best_rows.append((gid, gcount, 0, 0, 0, 0, 0, 0, "missed"))
            continue

        best = None
        for pid, inter in candidates:
            pcount = pred_count_dict.get(pid, int((pred == pid).sum()))
            union = gcount + pcount - inter
            iou = inter / union if union else 0.0
            dice = (2 * inter) / (gcount + pcount) if (gcount + pcount) else 0.0
            if best is None or iou > best["iou"]:
                best = {
                    "pid": pid,
                    "inter": inter,
                    "pcount": pcount,
                    "iou": iou,
                    "dice": dice,
                }

        status = "matched"
        for thr in matched_gt_at:
            if best["iou"] >= thr:
                matched_gt_at[thr] += 1

        gt_best_rows.append((
            gid,
            gcount,
            best["pid"],
            best["pcount"],
            best["inter"],
            best["iou"],
            best["dice"],
            len(candidates),
            status,
        ))

    write_csv(
        out_dir / "gt_best_match.csv",
        [
            "gt_id", "gt_voxels", "best_pred_id", "pred_voxels",
            "intersection", "iou", "dice", "num_pred_touching_gt", "status"
        ],
        gt_best_rows,
    )

    # Best GT for each predicted label.
    pred_best_rows = []
    for pid in pred_ids:
        pid = int(pid)
        pcount = pred_count_dict[pid]
        candidates = by_pred.get(pid, [])

        if not candidates:
            pred_best_rows.append((pid, pcount, 0, 0, 0, 0, 0, 0, "background_or_unmatched"))
            continue

        best = None
        valid_touch_count = 0

        for gid, inter in candidates:
            if gid not in valid_gt_set:
                continue
            valid_touch_count += 1
            gcount = gt_count_dict.get(gid, int((gt == gid).sum()))
            union = gcount + pcount - inter
            iou = inter / union if union else 0.0
            dice = (2 * inter) / (gcount + pcount) if (gcount + pcount) else 0.0
            if best is None or iou > best["iou"]:
                best = {
                    "gid": gid,
                    "inter": inter,
                    "gcount": gcount,
                    "iou": iou,
                    "dice": dice,
                }

        if best is None:
            pred_best_rows.append((pid, pcount, 0, 0, 0, 0, 0, valid_touch_count, "no_valid_gt"))
        else:
            status = "merge_candidate" if valid_touch_count > 1 else "matched_candidate"
            pred_best_rows.append((
                pid,
                pcount,
                best["gid"],
                best["gcount"],
                best["inter"],
                best["iou"],
                best["dice"],
                valid_touch_count,
                status,
            ))

    write_csv(
        out_dir / "pred_best_match.csv",
        [
            "pred_id", "pred_voxels", "best_gt_id", "gt_voxels",
            "intersection", "iou", "dice", "num_gt_touching_pred", "status"
        ],
        pred_best_rows,
    )

    # Marker diagnostics, optional.
    marker_components = 0
    marker_on_gt = None
    marker_on_bg = None

    if args.markers:
        markers = read3d(args.markers)
        marker_ids, marker_counts = label_counts(markers)
        marker_components = int(len(marker_ids))

        marker_mask = markers > 0
        marker_on_gt = int(np.logical_and(marker_mask, gt > 0).sum())
        marker_on_bg = int(np.logical_and(marker_mask, gt == 0).sum())

        marker_rows = []
        for mid in marker_ids:
            mid = int(mid)
            m = markers == mid
            touched = gt[m]
            touched = touched[touched > 0]
            if touched.size == 0:
                marker_rows.append((mid, int(m.sum()), 0, 0, "background"))
            else:
                ids, counts = np.unique(touched, return_counts=True)
                j = int(np.argmax(counts))
                marker_rows.append((mid, int(m.sum()), int(ids[j]), int(counts[j]), "inside_gt"))

        write_csv(
            out_dir / "marker_gt_overlap.csv",
            ["marker_id", "marker_voxels", "dominant_gt_id", "overlap_voxels", "status"],
            marker_rows,
        )

    thresholds = [float(x) for x in args.iou_thresholds.split(",")]

    summary_rows = [
        ("gt_total_ids_touching_crop", int(len(gt_ids_all))),
        ("gt_valid_ids_after_filter", int(len(valid_gt))),
        ("gt_border_ids_ignored_if_enabled", int(len(border_ids)) if args.ignore_gt_border else 0),
        ("gt_min_area_filter", int(args.min_gt_area)),
        ("ignore_gt_border", int(args.ignore_gt_border)),
        ("pred_label_count", int(len(pred_ids))),
        ("foreground_gt_voxels", gt_fg_count),
        ("foreground_pred_voxels", pred_fg_count),
        ("foreground_intersection", inter_fg),
        ("foreground_iou", fg_iou),
        ("foreground_dice", fg_dice),
        ("foreground_precision", fg_precision),
        ("foreground_recall", fg_recall),
        ("marker_components", marker_components),
    ]

    if marker_on_gt is not None:
        summary_rows.extend([
            ("marker_voxels_on_gt", marker_on_gt),
            ("marker_voxels_on_background", marker_on_bg),
        ])

    for thr in thresholds:
        matched = sum(1 for row in gt_best_rows if isinstance(row[5], float) and row[5] >= thr)
        recall = matched / len(valid_gt) if len(valid_gt) else 0.0
        summary_rows.append((f"gt_object_recall_iou_{thr}", recall))
        summary_rows.append((f"gt_matched_count_iou_{thr}", matched))

    write_csv(out_dir / "summary.csv", ["metric", "value"], summary_rows)

    print("")
    print("======================================")
    print("[eval] SUMMARY")
    print("======================================")
    for k, v in summary_rows:
        print(f"{k}: {v}")

    print("")
    print("[eval] wrote:")
    print(out_dir / "summary.csv")
    print(out_dir / "gt_best_match.csv")
    print(out_dir / "pred_best_match.csv")
    print(out_dir / "pair_intersections.csv")
    if args.markers:
        print(out_dir / "marker_gt_overlap.csv")


if __name__ == "__main__":
    main()
