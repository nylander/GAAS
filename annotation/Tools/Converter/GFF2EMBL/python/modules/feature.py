#!/usr/bin/env python2.7

import sys
from location import Location
from qualifiers import *

#
# required for next project:
#
#   NcRNA -> ncRNA
#   Non_canonical_five_prime_splice_site
#   Non_canonical_three_prime_splice_site
#   Rfam
#   Transcript -?-> does 'prim_transcript' work?


class Feature(object):
    """
    Super-class for the various feature types defined in:
    http://www.insdc.org/files/feature_table.html.
    """
    
    def __init__(self, feature, translate = {}, disregard = []):
        self.type = None
        self.location = Location(feature.location)
        self.translate = translate
        self.disregard = disregard
        self._load_default_translations()
        self.mandatory_qualifiers = {}
        self.optional_qualifiers = {}
    
    def __repr__(self):
        output = "\nFT   %s %s" % ("{:15}".format(self.type), self.location)
        for qualifier, value in self.mandatory_qualifiers.iteritems():
            if value:
                output += str(value)
        for qualifier, value in self.optional_qualifiers.iteritems():
            if value:
                output += str(value)
        return output
    
    def _load_default_translations(self):
        if "Dbxref" not in self.translate:
            self.translate["Dbxref"] = "db_xref"
        if "description" not in self.translate:
            self.translate["description"] = "note"
    
    @staticmethod
    def parse(feature, output = None):
        if output == None:
            output = []
        if feature.type == "gene":
            output += [GeneFeature(feature)]
        else:
            sys.stderr.write( "Unknown feature type '%s'" % feature.type + "\n" )
        #for sub_feature in feature.sub_features:
        #    output += self.parse(sub_feature)
        return output
    
    def add_qualifier(self, qualifier, value):
        if qualifier in self.disregard:
            return
        if qualifier in self.translate:
            qualifier = self.translate[qualifier]
        
        # handle mandatory qualifiers
        if qualifier not in self.mandatory_qualifiers:
            for i,q in enumerate([q.lower() for q in self.mandatory_qualifiers]):
                if q == qualifier.lower():
                    qualifier = self.mandatory_qualifiers.keys()[i]
                    break
        
        # handle optional qualifiers
        if qualifier not in self.optional_qualifiers:
            for i,q in enumerate([q.lower() for q in self.optional_qualifiers]):
                if q == qualifier.lower():
                    qualifier = self.optional_qualifiers.keys()[i]
                    break
        
        qualifier_type=qualifier[0].upper() + qualifier[1:] + "Qualifier"
        if qualifier in self.optional_qualifiers or qualifier in self.mandatory_qualifiers:
            
            try:
                if qualifier in self.optional_qualifiers:
                    if self.optional_qualifiers[qualifier] == None:
                        if type(value) == type(""):
                            self.optional_qualifiers[qualifier] = eval("%s('%s')" % (qualifier_type, value))
                        else:
                            self.optional_qualifiers[qualifier] = eval("%s(%s)" % (qualifier_type, value))
                    #else:
                    #    self.optional_qualifiers[qualifier].add(value)
                elif qualifier in self.mandatory_qualifiers:
                    if self.mandatory_qualifiers[qualifier] == None:
                        if type(value) == type(""):
                            self.mandatory_qualifiers[qualifier] = eval("%s('%s')" % (qualifier_type, value))
                        else:
                            self.mandatory_qualifiers[qualifier] = eval("%s(%s)" % (qualifier_type, value))
                    #else:
                    #    self.mandatory_qualifiers[qualifier].add(value)
            except Exception as e:
                sys.stderr.write( str(e) )
                import traceback
                traceback.print_exc(limit=5)
                
                sys.stderr.write( "Unknown qualifier '%s' for feature '%s'. " % (qualifier, type(qualifier)) + "\n" )
        else:
            sys.stderr.write( "unknown qualifier '%s' with value '%s' in %s" % (qualifier, value, type(self)) + "\n" )

# class Assembly_gapFeature
# class C_region
class CDSFeature( Feature ):

    definition="""coding sequence; sequence of nucleotides that
                  corresponds with the sequence of amino acids in a
                  protein (location includes stop codon);
                  feature includes amino acid conceptual translation."""

    comment = """/codon_start has valid value of 1 or 2 or 3, indicating
                 the offset at which the first complete codon of a coding
                 feature can be found, relative to the first base of
                 that feature;
                 /transl_table defines the genetic code table used if
                 other than the universal genetic code table;
                 genetic code exceptions outside the range of the specified
                 tables is reported in /transl_except qualifier;
                 /protein_id consists of a stable ID portion (3+5 format
                 with 3 position letters and 5 numbers) plus a version
                 number after the decimal point; when the protein
                 sequence encoded by the CDS changes, only the version
                 number of the /protein_id value is incremented; the
                 stable part of the /protein_id remains unchanged and as
                 a result will permanently be associated with a given
                 protein;"""

    def __init__(self, feature = None, translate = {'phase':'codon_start', "Name":"standard_name", "ID":"standard_name", "gene_name":"gene", "gene_type":"function"}, disregard = ["source", "Parent"]):
        super(CDSFeature, self).__init__(feature, translate, disregard)
        self.type = "CDS"

        self.optional_qualifiers = {'allele':None,
                           'artificial_location':None,
                           'citation':None,
                           'codon_start':None,
                           'db_xref':None,
                           'EC_number':None,
                           'exception':None,
                           'experiment':None,
                           'function':None,
                           'gene':None,
                           'gene_synonym':None,
                           'inference':None,
                           'locus_tag':None,
                           'map':None,
                           'note':None,
                           'number':None,
                           'old_locus_tag':None,
                           'operon':None,
                           'product':None,
                           'protein_id':None,
                           'pseudo':None,
                           'pseudogene':None,
                           'ribosomal_slippage':None,
                           'standard_name':None,
                           'translation':None,
                           'transl_except':None,
                           'transl_table':None,
                           'trans_splicing':None,
                           }

        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)

# class Centromere
# class D-loop
# class D_segment
class ExonFeature( Feature ):

    definition="""region of genome that codes for portion of spliced mRNA,
                  rRNA and tRNA; may contain 5'UTR, all CDSs and 3' UTR;"""

    def __init__(self, feature = None, translate = {"Name":"standard_name", "ID":"standard_name", "gene_name":"gene", "gene_type":"function"}, disregard = ['phase', "source", "Parent"]):
        super(ExonFeature, self).__init__(feature, translate, disregard)
        self.type = "exon"

        self.optional_qualifiers = {'allele':None,
                           'citation':None,
                           'db_xref':None,
                           'EC_number':None,
                           'experiment':None,
                           'function':None,
                           'gene':None,
                           'gene_synonym':None,
                           'inference':None,
                           'locus_tag':None,
                           'map':None,
                           'note':None,
                           'number':None,
                           'old_locus_tag':None,
                           'product':None,
                           'pseudo':None,
                           'pseudogene':None,
                           'standard_name':None,
                           'trans_splicing':None,
                           }

        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)

class GapFeature( Feature ):

    definition="""the location span of the gap feature for an unknown 
                  gap is 100 bp, with the 100 bp indicated as 100 "n"'s in 
                  the sequence.  Where estimated length is indicated by 
                  an integer, this is indicated by the same number of 
                  "n"'s in the sequence. 
                  No upper or lower limit is set on the size of the gap."""

    def __init__(self, feature = None, translate = {"Name":"standard_name", "ID":"standard_name"}, disregard = ['phase', "source", "Parent"]):
        super(GapFeature, self).__init__(feature, translate, disregard)
        self.type = "gap"
        
        self.mandatory_qualifiers = {'estimated_length':None,}
        
        self.optional_qualifiers = {'experiment':None,
                                    'inference':None,
                                    'locus_tag':None,
                                    'map':None,
                                    'note':None,
                                    }
        
        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)

class GeneFeature( Feature ):
    
    definition="""messenger RNA; includes 5'untranslated region (5'UTR),
                  coding sequences (CDS, exon) and 3'untranslated region
                  (3'UTR);"""
    
    comment = """a non-protein-coding gene, other than ribosomal RNA and
                 transfer RNA, the functional molecule of which is the RNA
                 transcript;"""
    
    def __init__(self, feature = None, translate = {"Name":"standard_name", "ID":"standard_name"}, disregard = ["orthology", "source", "phase", "eC_number"]):
        super(GeneFeature, self).__init__(feature, translate, disregard)
        self.type = "gene"
        self.optional_qualifiers = {'allele':None,
                           'artificial_location':None,
                           'citation':None,
                           'db_xref':None,
                           'experiment':None,
                           'function':None,
                           'gene':None,
                           'gene_synonym':None,
                           'inference':None,
                           'locus_tag':None,
                           'map':None,
                           'note':None,
                           'old_locus_tag':None,
                           'operon':None,
                           'product':None,
                           'pseudo':None,
                           'pseudogene':None,
                           'standard_name':None,
                           'trans_splicing':None,
                           }
        
        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)

# class IDNA
# class Intron
# class J_segment
# class LTR
# class Mat_peptide
# class Misc_binding
# class Misc_difference
# class Misc_feature
# class Misc_recomb
class Misc_RNAFeature( Feature ):

    definition="""any transcript or RNA product that cannot be defined by
                  other RNA keys (prim_transcript, precursor_RNA, mRNA,
                  5'UTR, 3'UTR, exon, CDS, sig_peptide, transit_peptide,
                  mat_peptide, intron, polyA_site, ncRNA, rRNA and tRNA);"""

    def __init__(self, feature = None, translate={"Name":"standard_name", "ID":"standard_name"}, disregard = ['phase', "source", "Parent"]):
        super(Misc_RNAFeature, self).__init__(feature, translate, disregard)
        self.type = "misc_RNA"
        self.optional_qualifiers = {'allele':None,
                           'citation':None,
                           'db_xref':None,
                           'experiment':None,
                           'function':None,
                           'gene':None,
                           'gene_synonym':None,
                           'inference':None,
                           'locus_tag':None,
                           'map':None,
                           'note':None,
                           'old_locus_tag':None,
                           'operon':None,
                           'product':None,
                           'pseudo':None,
                           'pseudogene':None,
                           'standard_name':None,
                           'trans_splitcing':None,
                           }

        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)

# class Misc_structure
# class Mobile_element
# class Modified_base
class MRNAFeature( Feature ):

    definition="""region of biological interest identified as a gene
                  and for which a name has been assigned"""

    def __init__(self, feature = None, translate={"Name":"standard_name", "ID":"standard_name"}, disregard = ['phase', "source", "eC_number", "Parent", "parent_name"]):
        super(MRNAFeature, self).__init__(feature, translate, disregard)
        self.type = "mRNA"
        self.optional_qualifiers = {'allele':None,
                           'citation':None,
                           'db_xref':None,
                           'experiment':None,
                           'function':None,
                           'gene':None,
                           'gene_synonym':None,
                           'inference':None,
                           'locus_tag':None,
                           'map':None,
                           'note':None,
                           'old_locus_tag':None,
                           'operon':None,
                           'product':None,
                           'pseudo':None,
                           'pseudogene':None,
                           'phenotype':None,
                           'standard_name':None,
                           'trans_splicing':None,
                           }

        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)

class NcRNAFeature( Feature ):

    definition="""a non-protein-coding gene, other than ribosomal RNA and
                  transfer RNA, the functional molecule of which is the RNA
                  transcript;"""
    
    comment="""the ncRNA feature is not used for ribosomal and transfer
               RNA annotation, for which the rRNA and tRNA feature keys
               should be used, respectively;"""
    
    def __init__(self, feature = None, translate={"Name":"standard_name", "ID":"standard_name"}, disregard = ['phase', "source", "eC_number", "Parent"]):
        super(NcRNAFeature, self).__init__(feature, translate, disregard)
        self.type = "mRNA"
        self.mandatory_qualifiers = {'ncRNA_class':None,
                                    }
        self.optional_qualifiers = {'allele':None,
                           'citation':None,
                           'db_xref':None,
                           'experiment':None,
                           'function':None,
                           'gene':None,
                           'gene_synonym':None,
                           'inference':None,
                           'locus_tag':None,
                           'map':None,
                           'note':None,
                           'old_locus_tag':None,
                           'operon':None,
                           'product':None,
                           'pseudo':None,
                           'pseudogene':None,
                           'standard_name':None,
                           'trans_splicing':None,
                           }

        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)



# class N_region
# class Old_sequence
# class Operon
# class OriT
# class PolyA_site
# class Precursor_RNA
# class Prim_transcript
# class Primer_bind
# class Protein_bind
# class Regulatory
# class Repeat_region
# class Rep_origin
class RRNAFeature( Feature ):

    definition="""mature ribosomal RNA; RNA component of the
                  ribonucleoprotein particle (ribosome) which assembles
                  amino acids into proteins"""

    def __init__(self, feature = None, translate={"Name":"standard_name", "ID":"standard_name"}, disregard = ['phase', "source", "Parent"]):
        super(RRNAFeature, self).__init__(feature, translate, disregard)
        self.type = "rRNA"
        self.optional_qualifiers = {'allele':None,
                           'citation':None,
                           'db_xref':None,
                           'experiment':None,
                           'function':None,
                           'gene':None,
                           'gene_synonym':None,
                           'inference':None,
                           'locus_tag':None,
                           'map':None,
                           'note':None,
                           'old_locus_tag':None,
                           'operon':None,
                           'product':None,
                           'pseudo':None,
                           'pseudogene':None,
                           'standard_name':None,
                           }

        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)

# class S_region
class Sig_peptideFeature( Feature ):

    definition="""signal peptide coding sequence; coding sequence for an
                  N-terminal domain of a secreted protein; this domain is
                  involved in attaching nascent polypeptide to the
                  membrane leader sequence;"""

    def __init__(self, feature = None, translate={"Name":"standard_name", "ID":"standard_name"}, disregard = ['phase', "source", "Parent"]):
        super(Sig_peptideFeature, self).__init__(feature, translate, disregard)
        self.type = "sig_peptide"

        self.optional_qualifiers = {'allele':None,
                           'citation':None,
                           'db_xref':None,
                           'experiment':None,
                           'function':None,
                           'gene':None,
                           'gene_synonym':None,
                           'inference':None,
                           'locus_tag':None,
                           'map':None,
                           'note':None,
                           'old_locus_tag':None,
                           'product':None,
                           'pseudo':None,
                           'pseudogene':None,
                           'standard_name':None,
                           }

        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)

class SourceFeature( Feature ):

    definition="""identifies the biological source of the specified span of
                  the sequence; this key is mandatory; more than one source
                  key per sequence is allowed; every entry/record will have, as a
                  minimum, either a single source key spanning the entire
                  sequence or multiple source keys, which together, span the
                  entire sequence."""

    def __init__(self, feature = None, translate={"Name":"standard_name", "ID":"standard_name"}, disregard = ['phase', "source"]):
        super(SourceFeature, self).__init__(feature, translate, disregard)
        self.type = "source"
        
        self.mandatory_qualifiers = {'organism':None,
                                     'mol_type':None,}
        
        self.optional_qualifiers = {'altitude':None,
                           'bio_material':None,
                           'cell_line':None,
                           'cell_type':None,
                           'chromosome':None,
                           'citation':None,
                           'clone':None,
                           'clone_lib':None,
                           'collected_by':None,
                           'collection_date':None,
                           'country':None,
                           'cultivar':None,
                           'culture_collection':None,
                           'db_xref':None,
                           'dev_stage':None,
                           'ecotype':None,
                           'environmental_sample':None,
                           'focus':None,
                           'germline':None,
                           'haplogroup':None,
                           'haplotype':None,
                           'host':None,
                           'identified_by':None,
                           'isolate':None,
                           'isolation_source':None,
                           'lab_host':None,
                           'lat_lon':None,
                           'macronuclear':None,
                           'map':None,
                           'mating_type':None,
                           'note':None,
                           'organelle':None,
                           'PCR_primers':None,
                           'plasmid':None,
                           'pop_variant':None,
                           'proviral':None,
                           'rearranged':None,
                           'segment':None,
                           'serotype':None,
                           'serovar':None,
                           'sex':None,
                           'specimen_voucher':None,
                           'strain':None,
                           'sub_clone':None,
                           'sub_species':None,
                           'sub_strain':None,
                           'tissue_lib':None,
                           'tissue_type':None,
                           'transgenic':None,
                           'type_material':None,
                           'variety':None,
                           }

        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)


# class Stem_loop
# class STS
# class Telomere
class TmRNAFeature( Feature ):

    definition="""transfer messenger RNA; tmRNA acts as a tRNA first,
                  and then as an mRNA that encodes a peptide tag; the
                  ribosome translates this mRNA region of tmRNA and attaches
                  the encoded peptide tag to the C-terminus of the
                  unfinished protein; this attached tag targets the protein for
                  destruction or proteolysis;"""

    def __init__(self, feature = None, translate={"Name":"standard_name", "ID":"standard_name"}, disregard = ['phase', "source", "eC_number", "Parent"]):
        super(TmRNAFeature, self).__init__(feature, translate, disregard)
        self.type = "tmRNA"
        self.optional_qualifiers = {'allele':None,
                           'citation':None,
                           'db_xref':None,
                           'experiment':None,
                           'function':None,
                           'gene':None,
                           'gene_synonym':None,
                           'inference':None,
                           'locus_tag':None,
                           'map':None,
                           'note':None,
                           'old_locus_tag':None,
                           'product':None,
                           'pseudo':None,
                           'pseudogene':None,
                           'standard_name':None,
                           'tag_peptide':None,
                           }

        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)

# class Transit_peptide
class TRNAFeature( Feature ):

    definition="""mature transfer RNA, a small RNA molecule (75-85 bases
                  long) that mediates the translation of a nucleic acid
                  sequence into an amino acid sequence;"""

    def __init__(self, feature = None, translate={"Name":"standard_name", "ID":"standard_name"}, disregard = ['phase', "source", "Parent"]):
        super(TRNAFeature, self).__init__(feature, translate, disregard)
        self.type = "tRNA"
        self.optional_qualifiers = {'allele':None,
                           'anticodon':None,
                           'citation':None,
                           'db_xref':None,
                           'experiment':None,
                           'function':None,
                           'gene':None,
                           'gene_synonym':None,
                           'inference':None,
                           'locus_tag':None,
                           'map':None,
                           'note':None,
                           'old_locus_tag':None,
                           'operon':None,
                           'product':None,
                           'pseudo':None,
                           'pseudogene':None,
                           'standard_name':None,
                           'trans_splicing':None,
                           }

        if feature:
            for qualifier, value in feature.qualifiers.iteritems():
                self.add_qualifier(qualifier, value)

# class Unsure
# class V_region
# class V_segment
# class Variation
# class 3'UTR
# class 5'UTR
#
# /!\ Specific features (UTR and CDS) can be "multiple line feature". It means, they are described through several lines if contain several exon.
# Depending of the tool used to produce the annotation, Some features may have several parents (like an exon share by multiple mRNA). It allows to factorise the information and condense the gff file.
# NBIS annotation service tool always "expand" them. Here the implementation don't manage that case. We report just an error that has to be fixed outside this code.
# The phylosophy of that code is expecting to have sequential describtion of feature ordered (i.e: all cds parts must be described in a row).
def parse_gff_feature(feature):
    
    features = []
    feature_type=feature.type[0].upper() + feature.type[1:] + "Feature"
        
    ### First check if a feature has a mulitple parents. In this case it through an error message.

    try:
        features += [eval("%s( feature )" % feature_type)] #### Where the different methods are called ####
    except Exception as e:
        sys.stderr.write( "WARNING parse_gff_feature: Unknown feature type '%s'\n" % feature_type )
    
    for sub_feature in feature.sub_features:
        features += parse_gff_feature(sub_feature)
    
    return features

if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser( description = __doc__ )
    parser.add_argument("gff_file", help="Annotation file in GFF3 format")
    args = parser.parse_args()
    
    try:
        from BCBio import GFF
        for record in GFF.parse( args.gff_file ):
            break
        
        for feature in record.features:
            for f in parse_gff_feature(feature):
                sys.stdout.write( f + "\n" )
            break
        
    except Exception as e:
        import traceback
        traceback.print_exc(limit=5)
        sys.stderr.write( str(e) )
        