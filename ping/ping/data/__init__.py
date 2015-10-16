from .base import filter_data, merge_by_key
from .default import PINGData
from .destrieux import DestrieuxData

prefixes = dict(desikan=dict(area='MRI_cort_area.ctx.', thickness='MRI_cort_thick.ctx.'),
                destrieux=dict(area='Destrieux_area.', thickness='Destrieux_thickness.'))
