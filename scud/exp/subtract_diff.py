from __future__ import division
from iotbx import mtz
from cctbx import miller
from scitbx.array_family import flex
import os,sys
from MTZ import MTZClass
import Plot

def run():
    '''
    read .mtz's and subtract everyting from the first one (in sys.argv)
    '''
    print("subtract_diff")
    mtz_names = sys.argv[1:]
    print(mtz_names)
    ### Generate descriptive output file name
    out_name = mtz_names[0][:-4]
    for fname in mtz_names[1:]:
        out_name+='-%s'%fname[:-4]
    out_name+='.mtz'
    ### First file is the starting file, from this, other arrays will be subtracted
    starting_mtz_object = mtz.object(mtz_names[0])
    starting_mtz_object.show_summary()
    ### Needed later for saving as .mtz file
    template = starting_mtz_object.as_miller_arrays()[2]
    ### Extrac data
    starting_miller_array_data = template.data()
    ### Loop over other .mtz files and subtract from starting data
    for mtz_file in mtz_names[1:]:
        temp_mtz_object = mtz.object(mtz_file)
        temp_mtz_object.show_summary()
        temp_miller_array_data = temp_mtz_object.as_miller_arrays()[2].data()
        starting_miller_array_data -= temp_miller_array_data

    ### Generate miller array containing new data
    final_miller = template.customized_copy(data=starting_miller_array_data)
    ### Save new .mtz
    mtz_dataset = final_miller.as_mtz_dataset(column_root_label='SUB',column_types='J')
    mtz_dataset.mtz_object().write(out_name)

if __name__ == run():
    run()
