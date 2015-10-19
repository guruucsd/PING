from .base import filter_data, merge_by_key
from .default import PINGData
from .destrieux import DestrieuxData

prefixes = dict(desikan=dict(area='MRI_cort_area.ctx.', thickness='MRI_cort_thick.ctx.'),
                destrieux=dict(area='Destrieux_area.', thickness='Destrieux_thickness.'))


def prefix2measure(prefix):
    for measure, p in prefixes[prefix2atlas(prefix)].items():
        if p == prefix:
            return measure
    return None


def prefix2atlas(prefix):
    for atlas, measures in prefixes.items():
        for measure, p in measures.items():
            if p == prefix:
                return atlas
    return None
