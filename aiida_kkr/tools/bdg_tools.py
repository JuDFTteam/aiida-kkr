"""
Helper tools that deal with the anomalous density of the BdG formalism in KKR
"""

from aiida.engine import calcfunction
from aiida.orm import FolderData


@calcfunction
def get_anomalous_density_data(retrieved, rename_files=None):
    """
    Extract anomalous density files from a retrieved folder of a KkrCalculation
    and copy into a new FolderData. This FolderData is then returned and can be
    used as the anomalous_density FolderData input to a KkrCalculation.

    :param retrieved: retrieved FolderData of a parent KKRhost BdG calculation where anomalous density files are stored (called `den_lm_ir.AAA.1.txt` where AAA is an atom index)

    :param rename_files: Optional Dict node where mappings of file names are defined. This is helpful if the atom index in the new structure is different from the original calculation (e.g. atom N of the original structure corresponds to atom M of the new structure, then the renaming dict should be {N: M}). Indices that are not found or that map to None are skipped and will not appear in the returned FolderData!

    :returns: anomalous_density FolderData which contains the (possibly renamed) anomalous density files
    """
    BdG_files = [i for i in retrieved.list_object_names() if 'den_lm_ir' in i]

    anomalous_density = FolderData()
    for fname in BdG_files:
        with retrieved.open(fname, 'r') as _fin:
            # default is to use the same name as in the input
            fname_out = fname

            # rename file, if natom is found in the renaming dict
            if rename_files is not None:
                # find atom index (should be integer)
                natom = int(fname.split('.')[1])
                # we pick the new atom index and recreate the filename
                natom_new = {int(k): v for k, v in rename_files.get_dict().items()}.get(natom)
                if natom_new is None:
                    # if no renaming is mapping is found
                    # we do not copy this anomalous density but leave it out
                    fname_out = None
                else:
                    # construct changed name
                    fname_out = 'den_lm_ir.%0.3i.1.txt' % natom_new

            if fname_out is not None:
                # copy input file to FolderData node, maybe with changed name
                anomalous_density.put_object_from_filelike(_fin, fname_out)

    return anomalous_density
