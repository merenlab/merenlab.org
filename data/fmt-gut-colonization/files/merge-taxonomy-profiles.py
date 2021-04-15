import sys
import argparse

from anvio.merger import MultipleRuns
from anvio.constants import levels_of_taxonomy
from anvio.utils import store_dict_as_TAB_delimited_file

m = MultipleRuns(argparse.Namespace())
m.input_profile_db_paths = sys.argv[1:-1]

m.populate_layer_additional_data_dict(missing_default_data_group_is_OK=True)

for level in levels_of_taxonomy:
    store_dict_as_TAB_delimited_file(m.layer_additional_data_dict[level],
                                     '%s_%s.txt' % (sys.argv[-1], level))
