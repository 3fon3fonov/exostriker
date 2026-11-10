"""Session migration for independent Hipparcos/Gaia likelihood selection."""


def ensure_dataset_selection(session):
    """Old sessions evaluated every loaded dataset regardless of legacy flags.

    Preserve that behavior once on migration; subsequent choices live in the
    existing dataset use flag and survive session serialization.
    """
    if getattr(session, '_absolute_astrometry_selection_version', 0) >= 1:
        return
    for name in ('ast_data_sets_hipp', 'ast_data_sets_gaia'):
        datasets = getattr(session, name, {})
        for dataset in datasets.values() if isinstance(datasets, dict) else datasets:
            if len(dataset):
                dataset[2] = True
    session._absolute_astrometry_selection_version = 1
