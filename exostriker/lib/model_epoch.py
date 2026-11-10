"""Resolve automatic model epochs without changing an explicitly chosen epoch."""
import numpy as np


def _first(sets):
    return min((float(np.min(row[0])) for row in sets.values() if len(row) and len(row[0])), default=None)


def first_data_epoch(session):
    """Automatic GUI epoch: RV observations take priority over photometry."""
    first_rv = _first(getattr(session, 'rv_data_sets', {}))
    return first_rv if first_rv is not None else _first(getattr(session, 'tra_data_sets', {}))


def resolve_epoch(session, include_rv=None, include_transits=None):
    types = getattr(session, 'type_fit', {})
    rv_sets = getattr(session, 'rv_data_sets', {})
    transit_sets = getattr(session, 'tra_data_sets', {})
    if include_rv is None:
        include_rv = types.get('RV', bool(rv_sets))
    if include_transits is None:
        include_transits = types.get('Transit', bool(transit_sets))
    epoch = float(session.epoch)
    # Transit import may seed the default before an RV dataset is added.
    automatic_transit = (hasattr(session, '_transit_default_epoch') and
                         epoch == session._transit_default_epoch)
    if epoch != 0. and not automatic_transit:
        return epoch
    first_rv = _first(rv_sets) if include_rv else None
    first_transit = _first(transit_sets) if include_transits else None
    if first_rv is not None and (epoch == 0. or automatic_transit):
        return first_rv
    if first_rv is None and first_transit is not None and (epoch == 0. or automatic_transit):
        return first_transit
    return epoch


def initialize_import_epoch(session):
    has_rv = any(len(row) for row in session.rv_data_sets.values())
    has_transit = any(len(row) for row in session.tra_data_sets.values())
    old = float(session.epoch)
    epoch = resolve_epoch(session, has_rv, has_transit)
    automatic = old == 0. or old == getattr(session, '_transit_default_epoch', None)
    session.epoch = epoch
    if automatic and has_rv:
        if hasattr(session, '_transit_default_epoch'):
            del session._transit_default_epoch
    elif automatic and has_transit:
        session._transit_default_epoch = epoch
