from MarkLib import PeptideList, Descriptor
from sys import argv
from math import log
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import percentileofscore
from os import path

inputfile = str(argv[1])
preFDR = 1
FDR = 1
RT_type = 'biolccc'
min_charge = 1
max_charge = 100
valid_proteins_input = '/home/mark/Python/48prots_new/48prots_percolator.txt'


def main(inputfile, preFDR, FDR, RT_type, min_charge, max_charge, valid_proteins_input):
    print 'inputfile = %s' % (inputfile, )
    line = inputfile
    peptides = PeptideList()
    peptides.get_from_pepxmlfile(line, min_charge=min_charge, max_charge=max_charge, max_rank=1)

    print 'total number of PSMs = %d' % (len(peptides.peptideslist),)

    if valid_proteins_input:
        valid_proteins = [line.strip() for line in open(valid_proteins_input, 'r')]
    else:
        valid_proteins = []

    # peptides.filter_titin()

    true_prots = set()
    prots_dict = {}
    pepts_dict = {}
    for peptide in peptides.peptideslist:
        try:
            pepts_dict[peptide.sequence] += 1
        except:
            pepts_dict[peptide.sequence] = 1
        for protein in peptide.parentproteins:
            try:
                prots_dict[protein.dbname] += 1
            except:
                prots_dict[protein.dbname] = 1
            if protein.dbname not in valid_proteins or peptide.note == 'decoy':
                protein.note = 'W'
                if peptide.note2 != 'tr' or peptide.note == 'decoy':
                    peptide.note2 = 'wr'
            else:
                true_prots.add(protein.dbname)
                protein.note = 'Valid'
                peptide.note2 = 'tr'

    for peptide in peptides.peptideslist:
        peptide.peptscore2 = pepts_dict[peptide.sequence]
        for protein in peptide.parentproteins:
            if peptide.protscore2 < float(prots_dict[protein.dbname]):
                peptide.protscore2 = float(prots_dict[protein.dbname])

    copy_peptides = PeptideList()
    copy_peptides.peptideslist = list(peptides.peptideslist)
    copy_peptides.pepxml_type = peptides.pepxml_type
    copy_peptides.filter_evalue_new(FDR=FDR, useMP=False)


    def PSMs_info(peptides):
        true_prots = set()
        for peptide in peptides.peptideslist:
            for protein in peptide.parentproteins:
                if protein.dbname in valid_proteins and peptide.note != 'decoy':
                    true_prots.add(protein.dbname)
        print 'True PSMs= %s' % (len([1 for x in peptides.peptideslist if x.note2 == 'tr']), )
        print 'Wrong PSMs= %s' % (len([1 for x in peptides.peptideslist if x.note2 == 'wr']), )
        print 'True prots = %s\n' % (len(true_prots))

    print 'Default filtering:'
    PSMs_info(copy_peptides)

    if RT_type == 'biolccc':
        copy_peptides.get_RC()
        peptides.RC = copy_peptides.RC
        copy_peptides.calc_RT(RTtype=RT_type)
        print copy_peptides.get_calibrate_coeff(RTtype=RT_type)
        peptides.calc_RT(RTtype=RT_type)
    else:
        copy_peptides.filter_modifications()
        copy_peptides.calc_RT(RTtype=RT_type)
        calibrate_coeff = copy_peptides.get_calibrate_coeff(RTtype=RT_type)
        copy_peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)
        copy_peptides.filter_RT(RT_tolerance=3 * calibrate_coeff[3])
        copy_peptides.calc_RT(RTtype=RT_type)
        calibrate_coeff = copy_peptides.get_calibrate_coeff(RTtype=RT_type)
        print calibrate_coeff
        peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)
        copy_peptides.peptideslist = list(peptides.peptideslist)
        copy_peptides.filter_evalue_new(FDR=preFDR, useMP=False)
        copy_peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)

    descriptors = []
    descriptors.append(Descriptor(name='RT difference',
        formula='(peptide.RT_exp - peptide.RT_predicted)',
        group='A', binsize='auto'))
    descriptors.append(Descriptor(name='precursor mass difference',
        formula='peptide.mass_diff()', group='A', binsize='auto'))
    descriptors.append(Descriptor(name='missed cleavages',
        formula='peptide.num_missed_cleavages', group='A', binsize='1'))
    descriptors.append(Descriptor(name='charge states',
        formula='peptide.pcharge', group='A', binsize='1'))
    descriptors.append(Descriptor(name='potential modifications',
        formula="(peptide.modified_code.count('[') - peptide.modified_code.count('[160]'))",
        group='A', binsize='1'))
    descriptors.append(Descriptor(name='isotopes mass difference',
        formula='round(peptide.massdiff, 0)', group='A', binsize='1'))
    descriptors.append(Descriptor(name='PSMs per protein',
        formula="peptide.protscore2", group='B'))
    descriptors.append(Descriptor(name='PSM count',
        formula='peptide.peptscore2', group='B'))


    def plot_histograms(descriptors, peptides):
        fig = plt.figure(figsize=(16, 12))

        def find_optimal_xy(descriptors):
            x, y = 1, 1
            while x * y < len(descriptors) + 1:
                if x > y:
                    y += 1
                else:
                    x += 1
            return x, y
        ox, oy = find_optimal_xy(descriptors)
        copy_peptides = PeptideList()
        copy_peptides.peptideslist = list(peptides.peptideslist)
        copy_peptides.pepxml_type = peptides.pepxml_type
        threshold1, threshold2 = copy_peptides.filter_evalue_new(FDR=FDR, useMP=True)
        threshold1 = -log(threshold1)
        threshold2 = log(threshold2)
        for idx, descriptor in enumerate(descriptors):
            ax = fig.add_subplot(ox, oy, idx + 1)
            array_wrong = [eval(descriptor.formula)
                    for peptide in peptides.peptideslist if peptide.note2 == 'wr']
            array_valid = [eval(descriptor.formula)
                    for peptide in peptides.peptideslist if peptide.note2 == 'tr']
            if descriptor.group == 'B':
                array_wrong = np.log10(array_wrong)
                array_valid = np.log10(array_valid)
            binsize = float(descriptor.get_binsize(copy_peptides))
            lbin = min(np.append(array_wrong, array_valid))
            rbin = max(np.append(array_wrong, array_valid)) + 2 * binsize
            if rbin == lbin:
                rbin += binsize
            H1, _ = np.histogram(array_wrong, bins=np.arange(lbin, rbin, binsize))
            H2, _ = np.histogram(array_valid, bins=np.arange(lbin, rbin, binsize))
            H1 = H1.clip(1)
            H2 = H2.clip(1)
            if descriptor.group == 'B':
                H1 = np.log10(H1)
                H2 = np.log10(H2)
            ind = np.arange(lbin, rbin - 1.01 * binsize, binsize)
            width = binsize
            plt.bar(ind, H1, width, color='red', alpha=0.5)
            plt.bar(ind, H2, width, color='green', alpha=0.5)
            ax.set_xlabel(descriptor.name)

        zero_peptscore = log(min(
            p.peptscore for p in peptides.peptideslist if p.peptscore != 0)) - 1
        zero_evalue = -log(min(
            p.evalue for p in peptides.peptideslist if p.evalue != 0)) + 1
        PSMs_wrong = [[(-log(pept.evalue) if pept.evalue != 0 else zero_evalue),
            (log(pept.peptscore) if pept.peptscore != 0 else zero_peptscore)]
            for pept in peptides.peptideslist if pept.note2 == 'wr']
        PSMs_true = [[(-log(pept.evalue) if pept.evalue != 0 else zero_evalue),
            (log(pept.peptscore) if pept.peptscore != 0 else zero_peptscore)]
            for pept in peptides.peptideslist if pept.note2 == 'tr']

        print 'MP filtering:'
        PSMs_info(copy_peptides)
        print 'Without filtering:'
        PSMs_info(peptides)

        ax = fig.add_subplot(ox, oy, 9)
        ax.plot([x[0] for x in PSMs_wrong], [x[1] for x in PSMs_wrong], 'o',
                markersize=2, color='red')
        ax.plot([x[0] for x in PSMs_true], [x[1] for x in PSMs_true], 'o',
                markersize=2, color='blue')
        ax.axvline(threshold1, color='green')
        if threshold2:
            ax.axhline(threshold2, color='green')
        ax.set_ylim(min(x[1] for x in PSMs_true) - 1, max(x[1] for x in PSMs_true) + 1)
        fig.tight_layout()
        plt.show()
#       plt.savefig('%s/Results_new_%s.png' % (path.dirname(path.realpath(inputfile)), path.splitext(path.splitext(path.basename(inputfile))[0])[0]))


    def prepare_hist(descriptors, copy_peptides):
        for descriptor in descriptors:
            descriptor.array = descriptor.get_array(copy_peptides)
            if descriptor.group == 'A':
                binsize = float(descriptor.get_binsize(copy_peptides))
                lbin, rbin = min(descriptor.array), max(descriptor.array)
                if lbin == rbin:
                    rbin = lbin + binsize
                descriptor.hist = np.histogram(
                        descriptor.array, bins=np.arange(lbin, rbin, binsize))
        return descriptors


    def calc_peptscore(peptide, descriptors):
        for descriptor in descriptors:
            if descriptor.group == 'A':
                j = len(descriptor.hist[1]) - 2
                if j > -1:
                    while j >= -1:
                        if j == -1:
                            peptide.peptscore = 0
                            break
                        if descriptor.hist[1][j] <= eval(descriptor.formula):
                            peptide.peptscore *= float(descriptor.hist[0][j]) / sum(descriptor.hist[0])
                            break
                        j -= 1
            elif descriptor.group == 'B':
                peptide.peptscore *= percentileofscore(descriptor.array, eval(descriptor.formula)) / 100

    descriptors = prepare_hist(descriptors, copy_peptides)
    for peptide in peptides.peptideslist:
        calc_peptscore(peptide, descriptors)
    plot_histograms(descriptors, peptides)

if __name__ == '__main__':
    main(inputfile, preFDR, FDR, RT_type,
            min_charge, max_charge, valid_proteins_input)
