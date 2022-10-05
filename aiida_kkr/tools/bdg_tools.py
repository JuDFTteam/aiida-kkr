"""
Helper tools that deal with the anomalous density of the BdG formalism in KKR
"""

from aiida.engine import calcfunction


@calcfunction
def get_anomalous_density_data(retrieved):
    """
    Extract anomalous density files from a retrieved folder of a KkrCalculation
    and copy into a new FolderData. This FolderData is then returned and can be
    used as the anomalous_density FolderData input to a KkrCalculation.
    """
    BdG_files = [i for i in retrieved.list_object_names() if 'den_lm_ir' in i]

    anomalous_density = orm.FolderData()
    for fname in BdG_files:
        with retrieved.open(fname, 'r') as _fin:
            anomalous_density.put_object_from_filelike(_fin, fname)

    return anomalous_density
