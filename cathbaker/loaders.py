# core
import json
import logging
import os

# non-core
import requests
from cathpy.align import Align

# local
from cathbaker import models

LOG = logging.getLogger(__name__)

API_FUNFAM_BASE = 'https://www.cathdb.info/version/{cath_version}/api/rest/funfam'
API_FUNFAM_ALN = 'https://www.cathdb.info/version/{cath_version}/api/rest/superfamily/{superfamily_id}/funfam/{funfam_number}/files/seed_alignment'

def get_funfam_alignment(*, cath_version, superfamily_id, funfam_number, **kwargs):
    """
    Given a FunfamInfo
    """
    ff_url = API_FUNFAM_ALN.format(cath_version=cath_version, superfamily_id=superfamily_id, funfam_number=funfam_number)
    
    LOG.info("GET %s", ff_url)
    response = requests.get(ff_url)
    response.raise_for_status()

    aln = Align.from_stockholm(response.text)
    return aln

def load_all_funfam_info(cath_version, *, cache_file=None, nocache=False):
    """
    Loads Funfam info from API
    """

    funfam_url = API_FUNFAM_BASE.format(cath_version=cath_version)

    if not cache_file:
        cache_file = f'cath.{cath_version}.funfam.json'

    # write the API json data to a local cache
    if not os.path.isfile(cache_file) or nocache:

        LOG.info(f"Fetching Funfam data from '{funfam_url}' ...")
        response = requests.get(funfam_url)
        response.raise_for_status()
        LOG.info("   ... done")

        LOG.info(f"Writing to cache file '{cache_file}'")
        with open(cache_file, 'wt') as fh:
            fh.write(response.text)

    # read the json
    with open(cache_file) as fh:
        LOG.info(f"Loading cached Funfam data from '{cache_file}' ...")
        all_funfam_data = json.load(fh)

    funfams = []
    for ff_data in all_funfam_data['data']:
        ff = models.FunfamInfo(**ff_data, cath_version=cath_version)
        funfams.extend([ff])

    return funfams
