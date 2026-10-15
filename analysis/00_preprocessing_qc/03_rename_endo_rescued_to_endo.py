"""Rename ENDO_rescued_qc → ENDO in the merged atlas+endo h5ad's categorical
columns. Uses h5py for in-place edits — does NOT rewrite the X matrix.

Touches only the categorical metadata for these annotation columns:
  Celltype_markers
  subcluster_annot
  subcluster_annot_named
  subcluster_annot_labeled
  subcluster_annot_labeled_renamed
  subcluster_annot_labeled_renamed_neurons
  subcluster_annot_labeled_renamed_neurons2
  label_transfer_annot
  predictions_combined
  predictions_scanvi
  predictions_scanvi_grouped
  parent_cell_type
  parent_type

Leaves untouched: cell_origin, major_type (rescued cells stay 'VASCULAR'),
mic7_flag, raw_count, normalized, X, obsm, etc.
"""

import os
import h5py
import numpy as np

BASE = os.getcwd()  # Set to project root (containing h5ad/) before running
H5 = os.path.join(BASE, "h5ad", "13ad_3m_5k_atlas_plus_endo.h5ad")

COLS_TO_FIX = [
    "Celltype_markers",
    "subcluster_annot",
    "subcluster_annot_named",
    "subcluster_annot_labeled",
    "subcluster_annot_labeled_renamed",
    "subcluster_annot_labeled_renamed_neurons",
    "subcluster_annot_labeled_renamed_neurons2",
    "subcluster_annot_labeled_no_index",
    "subcluster_annot_labeled_no_index_with_pct",
    "label_transfer_annot",
    "predictions_combined",
    "predictions_scanvi",
    "predictions_scanvi_grouped",
    "parent_cell_type",
    "parent_type",
]

NEW_LABEL = "ENDO"
OLD_LABEL = "ENDO_rescued_qc"


def merge_categorical(grp, old, new):
    """Merge `old` category into `new` in an anndata-encoded categorical h5py group.

    Group structure: grp/categories (str array), grp/codes (int array).
    """
    cats_ds = grp["categories"]
    codes_ds = grp["codes"]
    cats = cats_ds[:].astype(str)
    codes = codes_ds[:].astype(np.int64)

    if old not in cats:
        return False, f"  {grp.name}: '{old}' not present (cats: {cats[:30].tolist()}...)"

    old_idx = int(np.where(cats == old)[0][0])

    if new not in cats:
        # Just rename the existing category in place
        cats[old_idx] = new
        del grp["categories"]
        grp.create_dataset("categories", data=np.array(cats, dtype=object),
                           dtype=h5py.string_dtype(encoding="utf-8"))
        return True, f"  {grp.name}: renamed '{old}' -> '{new}' (no merge)"

    new_idx = int(np.where(cats == new)[0][0])

    # Re-map old code -> new code
    codes[codes == old_idx] = new_idx

    # Remove the old category and shift codes greater than old_idx down by 1
    new_cats = np.delete(cats, old_idx)
    shift_mask = codes > old_idx
    codes[shift_mask] = codes[shift_mask] - 1

    # Write back
    del grp["categories"]
    grp.create_dataset(
        "categories",
        data=np.array(new_cats, dtype=object),
        dtype=h5py.string_dtype(encoding="utf-8"),
    )
    codes_ds[...] = codes.astype(codes_ds.dtype)
    return True, f"  {grp.name}: merged '{old}' (n_cells_relabeled={int(shift_mask.sum() + (codes == new_idx).sum() - shift_mask.sum())}) into '{new}'"


def main():
    print(f"Opening {H5} for in-place edits", flush=True)
    with h5py.File(H5, "r+") as f:
        obs = f["obs"]
        print(f"obs keys (first 30): {list(obs.keys())[:30]}", flush=True)
        for col in COLS_TO_FIX:
            if col not in obs:
                print(f"  [skip] '{col}' not in obs")
                continue
            grp = obs[col]
            if not isinstance(grp, h5py.Group):
                print(f"  [skip] '{col}' is not a categorical group")
                continue
            if "categories" not in grp or "codes" not in grp:
                print(f"  [skip] '{col}' lacks categories/codes")
                continue
            ok, msg = merge_categorical(grp, OLD_LABEL, NEW_LABEL)
            print(f"  [{'OK' if ok else '--'}] {msg}")
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
