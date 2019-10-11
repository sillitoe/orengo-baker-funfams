# core
import logging

LOG = logging.getLogger(__name__)

class FunfamInfo:
    """
    Stores information from Funfam API
    """

    def __init__(self, *, cath_version, name, funfam_number, superfamily_id, rep_id, rep_source_id, 
        seed_dops_score, inclusion_bitscore, inclusion_e_value,
        num_members_in_funfam, num_members_in_seed_aln, **kwargs):

        if seed_dops_score is None:
            seed_dops_score = 0
        if inclusion_bitscore is None:
            inclusion_bitscore = 0
        if inclusion_e_value is None:
            inclusion_e_value = 0
        if num_members_in_seed_aln is None:
            num_members_in_seed_aln = 0
        if num_members_in_funfam is None:
            num_members_in_funfam = 0

        self.cath_version = cath_version

        self.name = name
        self.funfam_number = int(funfam_number)
        self.superfamily_id = superfamily_id
        self.rep_id = rep_id
        self.rep_source_id = rep_source_id
        self.seed_dops_score = float(seed_dops_score)
        self.inclusion_bitscore = float(inclusion_bitscore)
        self.inclusion_e_value = float(inclusion_e_value)
        self.num_members_in_funfam = int(num_members_in_funfam)
        self.num_members_in_seed_aln = int(num_members_in_seed_aln)

    @property
    def topology_id(self):
        """Return the C.A.T code"""
        return '.'.join(self.superfamily_id.split('.')[:3])

    @property
    def arch_id(self):
        """Return the C.A code"""
        return '.'.join(self.superfamily_id.split('.')[:2])

    def __str__(self):
        return f'{self.superfamily_id}/FF/{str(self.funfam_number)}: "{self.name}" (DOPS: {str(self.seed_dops_score)})'
