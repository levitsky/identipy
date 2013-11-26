from MarkLib import PeptideList, Descriptor
from sys import argv
from math import log
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import percentileofscore
from os import path
import pickle
from pyteomics import mzml, fasta, auxiliary, mgf

inputfile = str(argv[1])
preFDR = 1
FDR = 1
RT_type = 'biolccc'
min_charge = 1
max_charge = 10

protsC = {}
protsL = {}
valid_proteins = []
valid_proteins_input = ''
if valid_proteins_input:
    with open(valid_proteins_input) as f:
        for line in f:
            lss = line.strip().split()
            if len(lss) > 1:
                dbname, conc = lss
                valid_proteins.append(dbname)
                protsC[dbname] = float(conc)
            else:
                dbname = line.strip().split()[0]
                valid_proteins.append(dbname)
#    valid_proteins = [line.strip() for line in open(valid_proteins_input)]
else:
    valid_proteins = []


def find_optimal_xy(descriptors):
    x, y = 1, 1
    while x * y < len(descriptors) + 1:
        if x > y:
            y += 1
        else:
            x += 1
    return x, y

def PSMs_info(peptides, tofile=False):
    proteins_groups = []
    prots = {}
    prots_peptides = {}
    peptides_added = set()
    true_prots = set()
    unique = set()
    todel = set()
    Total_prots = set()
    Total_prots_2 = set()
    for peptide in peptides.peptideslist:
        flag = 1
        if peptide.note2 == 'tr':
            unique.add(peptide.sequence)
        for group in proteins_groups:
            if any(protein.dbname in group for protein in peptide.parentproteins):
                for protein in peptide.parentproteins:
                    group.add(protein.dbname)
                flag = 0
                break
        if flag:
            proteins_groups.append(set(protein.dbname for protein in peptide.parentproteins))
        #Count PSMs:
        for protein in peptide.parentproteins:
            if protein.dbname in Total_prots:
                Total_prots_2.add(protein.dbname)
            else:
                Total_prots.add(protein.dbname)
            try:
                prots[protein.dbname]['PSMs'] += 1
                prots[protein.dbname]['sumI'] += peptide.sumI
            except:
                prots[protein.dbname] = {}
                prots[protein.dbname]['PSMs'] = 1
                prots[protein.dbname]['sumI'] = peptide.sumI
            if protein.dbname in valid_proteins and peptide.note != 'decoy':
                true_prots.add(protein.dbname)
        #Count Peptides:
        if peptide.sequence not in peptides_added:
            for protein in peptide.parentproteins:
                try:
                    prots[protein.dbname]['Peptides'] += 1
                except:
                    prots[protein.dbname]['Peptides'] = 1
        peptides_added.add(peptide.sequence)
    for group in proteins_groups:
        bestname = min(dbname for dbname in group
                if prots[dbname]['sumI'] == max(
                    prots[dn]['sumI'] for dn in group))
        for dbname in group:
            if dbname != bestname or prots[dbname]['PSMs'] < 2:
                todel.add(dbname)#del prots[dbname]
    if 'P02768' in prots.keys():
        print prots['P02768']['sumI']
    for dbname in todel:
        del prots[dbname]
    normK = 0
    for dbname in prots.keys():
        prots[dbname]['sumI'] /= protsL[dbname]
        normK += prots[dbname]['sumI']
    if 'P02768' in prots.keys():
        print prots['P02768']['sumI']
    for dbname in prots.keys():
        prots[dbname]['sumI'] /= normK
    if 'P02768' in prots.keys():
        print prots['P02768']['sumI']
    if tofile:
        output_proteins = open('%s/Results_new_%s_proteins.csv' % (
            path.dirname(path.realpath(inputfile)),
            path.splitext(path.splitext(path.basename(inputfile))[0])[0]), 'w')
        output_peptides = open('%s/Results_new_%s_eptides.csv' % (
            path.dirname(path.realpath(inputfile)),
            path.splitext(path.splitext(path.basename(inputfile))[0])[0]), 'w')
        output_peptides_detailed = open('%s/%s_peptides.csv' % (
            path.dirname(path.realpath(inputfile)),
            path.splitext(path.splitext(path.basename(inputfile))[0])[0]), 'w')
        output_peptides_detailed.write(
                'sequence\tmodified_sequence\te-value\tMPscore\tspectrum_title'
                '\tproteins\n')
        if protsC:
            output_proteins_valid = open(
                    '%s/Results_new_%s_proteins_valid.csv' % (
                        path.dirname(path.realpath(inputfile)),
                        path.splitext(
                            path.splitext(
                                path.basename(inputfile))[0])[0]), 'w')
            temp_data = []
        for k, v in prots.items():
            if protsC and k in valid_proteins:
                output_proteins_valid.write('%s,%s,%s,%s,%s\n' % (
                    k, v['PSMs'], v['Peptides'], v['sumI'], protsC[k]))
                temp_data.append([float(v['sumI']), protsC[k]])
            output_proteins.write('%s\t%s\t%s\t%s\t%s\n' % (
                k, v['PSMs'], v['Peptides'], v['sumI'],
                ('+' if k in valid_proteins else '-')))
        for peptide in peptides.peptideslist:
            if any(protein.dbname in
                    [k for k, v in prots.items() if int(v['PSMs']) > 1]
                    for protein in peptide.parentproteins):
#               if peptide.modified_code.count('['
#                      ) == peptide.modified_code.count('[160]'):
                output_peptides.write('%s\t%s\t%s\n' % (
                    peptide.sequence, peptide.evalue, peptide.RT_exp))
#                    print peptide.modified_code

        peptides_best = {}
        for peptide in peptides.peptideslist:
            if (peptide.sequence in peptides_best and
                    peptide.evalue < peptides_best[peptide.sequence]):
                peptides_best[peptide.sequence] = peptide.evalue
            else:
                peptides_best[peptide.sequence] = peptide.evalue
        for peptide in peptides.peptideslist:
            if peptide.evalue == peptides_best[peptide.sequence]:
                output_peptides_detailed.write('%s\t%s\t%s\t%s\t%s\t' % (
                    peptide.sequence, peptide.modified_code, peptide.evalue,
                    peptide.peptscore, peptide.spectrum))
                flag = 0
                for protein in peptide.parentproteins:
                    if flag == 0:
                        output_peptides_detailed.write('"%s"' % (protein.dbname))
                        flag = 1
                    else:
                        output_peptides_detailed.write(',"%s"' % (protein.dbname))
                output_peptides_detailed.write('\n')
        if protsC:
            temp_sum = sum([x[0] for x in temp_data])
            temp_data = [[x[0] / temp_sum, x[1]] for x in temp_data]
            print 'conc: ', auxiliary.linear_regression(
                    [x[0] for x in temp_data], [x[1] for x in temp_data])

        output_proteins.close()
    print 'unique = %d' % (len(unique))
    print 'True PSMs=', sum(1 for x in peptides.peptideslist if x.note2 == 'tr')
    print 'Wrong PSMs=', sum(1 for x in peptides.peptideslist if x.note2 == 'wr')
    print 'Prots =', len(prots)
    print 'Prots >= 2 =', sum(1 for v in prots.values() if v['PSMs'] > 1)
    print 'Total_prots >=2 =', len(Total_prots_2)
    print 'True Prots = {}\n'.format(len(true_prots))

def plot_histograms(descriptors, peptides):
    fig = plt.figure(figsize=(16, 12))
    ox, oy = find_optimal_xy(descriptors)
    copy_peptides = PeptideList()
    copy_peptides.peptideslist = list(peptides.peptideslist)
    copy_peptides.pepxml_type = peptides.pepxml_type
    print len(copy_peptides.peptideslist), '!!!!!'
    print len(peptides.peptideslist), '!!!!!'
    copy_peptides.filter_evalue_new(FDR=FDR, useMP=False)
    print len(copy_peptides.peptideslist), '!!!!!'
    for idx, descriptor in enumerate(descriptors):
        ax = fig.add_subplot(ox, oy, idx + 1)
        array_wrong = [eval(descriptor.formula)
                for peptide in peptides.peptideslist if peptide.note2 == 'wr']
                # and peptide.note == 'decoy']
        array_valid = [eval(descriptor.formula)
                for peptide in peptides.peptideslist if peptide.note2 == 'tr']
        if descriptor.group == 'B':
            array_wrong = np.log10(array_wrong)
            array_valid = np.log10(array_valid)
        binsize = float(descriptor.get_binsize(copy_peptides))
        if binsize < float(max(np.append(array_wrong, array_valid)) -
                min(np.append(array_wrong, array_valid))) / 400:
            binsize = float(max(np.append(array_wrong, array_valid)) -
                    min(np.append(array_wrong, array_valid))) / 400
        lbin = min(np.append(array_wrong, array_valid))
        rbin = max(np.append(array_wrong, array_valid)) + 1.5 * binsize
        if rbin == lbin:
            rbin += binsize
        if descriptor.name == 'fragment mass tolerance, Da':
            rbin = 0.05
        print descriptor.name
        print len(array_wrong), lbin, rbin, binsize
        H1, _ = np.histogram(array_wrong, bins=np.arange(lbin, rbin, binsize))
        H2, _ = np.histogram(array_valid, bins=np.arange(lbin, rbin, binsize))
        if descriptor.group == 'B':
            H3, _ = np.histogram([np.log10(eval(descriptor.formula))
                for peptide in copy_peptides.peptideslist],
                bins=np.arange(lbin, rbin, binsize))
        else:
            H3, _ = np.histogram([eval(descriptor.formula)
                for peptide in copy_peptides.peptideslist],
                bins=np.arange(lbin, rbin, binsize))
        if descriptor.name in ['precursor mass difference, ppm']:
            print '!!!!!!!!!!!!!!!!!! MEDIAN prec mass =', np.median(
                    [eval(descriptor.formula)
                        for peptide in copy_peptides.peptideslist])
        if descriptor.group == 'B':
            H1 = H1.clip(1)
            H2 = H2.clip(1)
            H3 = H3.clip(1)
            H1 = np.log10(H1)
            H2 = np.log10(H2)
            H3 = np.log10(H3)
        ind = np.arange(lbin, rbin - (1 + 1e-15) * binsize, binsize)
        width = binsize
        print descriptor.name, len(H1), len(ind)
        plt.bar(ind, H1, width, color='red', alpha=0.5)
        plt.bar(ind, H2, width, color='green', alpha=0.5)
        plt.bar(ind, H3, width, color='black', alpha=0.5)

#        if descriptor.name in ['precursor mass difference, ppm']:
#            H3, _ = np.histogram([eval(descriptor.formula)
#                for peptide in copy_peptides.peptideslist],
#                bins=np.arange(lbin, rbin, binsize))
#            plt.bar(ind, H3, width, color='black', alpha=0.5)
        if descriptor.name in ['missed cleavages', 'charge states',
                'potential modifications', 'isotopes mass difference, Da']:
            plt.xticks(range(1, 6))
        ax.set_xlabel(descriptor.name)
    return fig


def plot_MP(descriptors, peptides, fig):
    ox, oy = find_optimal_xy(descriptors)
    copy_peptides = PeptideList()
    copy_peptides.peptideslist = list(peptides.peptideslist)
    copy_peptides.pepxml_type = peptides.pepxml_type
    threshold1, threshold2 = copy_peptides.filter_evalue_new(FDR=FDR, useMP=True)
    threshold1 = -log(threshold1)
    try:
        threshold2 = log(threshold2)
    except:
        pass

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
    PSMs_info(copy_peptides, tofile=True)
    print 'Without filtering:'
    PSMs_info(peptides)

    ax = fig.add_subplot(ox, oy, ox * oy)
    ax.plot([x[0] for x in PSMs_wrong], [x[1] for x in PSMs_wrong],
            'o', markersize=2, color='red')
    ax.plot([x[0] for x in PSMs_true], [x[1] for x in PSMs_true],
            'o', markersize=2, color='blue')
    ax.axvline(threshold1, color='green')
    if threshold2:
        ax.axhline(threshold2, color='green')
    ax.set_ylim(min(x[1] for x in PSMs_true) - 1, max(x[1] for x in PSMs_true) + 1)
    fig.tight_layout()
#        plt.show()
    plt.savefig('%s/Results_new_%s.png' % (path.dirname(path.realpath(inputfile)),
        path.splitext(path.splitext(path.basename(inputfile))[0])[0]))

def prepare_hist(descriptors, copy_peptides):
    for descriptor in descriptors:
        descriptor.array = descriptor.get_array(copy_peptides)
        if descriptor.group == 'A':
            binsize = float(descriptor.get_binsize(copy_peptides))
            if binsize < float(max(descriptor.get_array(copy_peptides)) -
                    min(descriptor.get_array(copy_peptides))) / 400:
                    binsize = float(max(descriptor.get_array(copy_peptides)) -
                            min(descriptor.get_array(copy_peptides))) / 400

            lbin, rbin = min(descriptor.array), max(descriptor.array) + 1.5 * binsize
    #                if descriptor.name == 'RT difference, min':
#                    rbin = 100
            if lbin == rbin:
                rbin = lbin + binsize
            descriptor.hist = np.histogram(descriptor.array,
                    bins=np.arange(lbin, rbin + binsize, binsize))
    return descriptors


def calc_peptscore(peptide, descriptors, jk):
    for descriptor in descriptors:
        if peptide.peptscore != 0:
            flag = 1
        else:
            flag = 0
        if descriptor.group == 'A':
            if descriptor.hist[1][-1] < eval(descriptor.formula):
                peptide.peptscore = 0
            else:
                j = len(descriptor.hist[1]) - 2
                if j > -1:
                    while j >= -1:
                        if j == -1:
                            peptide.peptscore = 0
                            break
                        if descriptor.hist[1][j] <= eval(descriptor.formula):
                            peptide.peptscore *= float(
                                    descriptor.hist[0][j]) / sum(
                                            descriptor.hist[0])
                            break
                        j -= 1
        elif descriptor.group == 'B':
#                peptide.peptscore *= max(percentileofscore(descriptor.array,
#        eval(descriptor.formula)), percentileofscore(descriptor.array,
#        min(descriptor.array))) / 100
            peptide.peptscore *= percentileofscore(
                    descriptor.array, eval(descriptor.formula)) / 100

        if flag and not peptide.peptscore:
            try:
                jk[descriptor.name] += 1
            except:
                jk[descriptor.name] = 1

    return jk



def main(inputfile, preFDR, FDR, RT_type, min_charge, max_charge):
    print 'inputfile = %s' % (inputfile, )
    line = inputfile
    peptides = PeptideList()
    peptides.get_from_pepxmlfile(line, min_charge=min_charge,
            max_charge=max_charge, max_rank=1)
    print len(peptides.peptideslist), 'with modifications'
#    peptides.filter_modifications_test()
    print len(peptides.peptideslist), 'without modifications'
    """
    for peptide in peptides.peptideslist:
        if peptide.mass_diff() > 10 or peptide.mass_diff() < -10:
            print peptide.spectrum
            print peptide.sequence
            print peptide.mass_diff()
            print '\n'
    """
    if len(argv) > 2:
        fastafile = fasta.read(argv[2])
        for x in fastafile:
            if not any(x[0].startswith(tag)
                    for tag in ['sp', 'tr', 'DECOY_sp', 'DECOY_tr']):
                if any(tag in x[0] for tag in ['SWISS-PROT:', 'TREMBL:']):
                    dbname = x[0].split(' ')[0]
                else:
                    dbname = x[0]
                protsL[dbname] = len(x[1])
            else:
                protsL[x[0].split('|')[1]] = len(x[1])
    Fragment_intensities = {}
    if len(argv) == 3:
        argv.append(argv[1].split('.pep.xml')[0]+'.t.xml')
        try:
            argv.append(argv[1].split('.pep.xml')[0]+'.mgf')
        except:
            pass
    if len(argv) > 3:
        if argv[3].split('.')[-2:] == ['t', 'xml']:
            for x in open(argv[3], 'r').readlines():
                if x.startswith('<group id='):
                    Fragment_intensities[
                            int(x.split('<group id="')[1].split('"')[0])
                            ] = 10**float(x.split('sumI="')[1].split('"')[0])
        else:
            spectra = [x for x in mzml.read(argv[3]) if x['ms level'] == 1]
            for psm in peptides.peptideslist:
                j = len(spectra) - 1
                while spectra[j]['scanList']['scan'][0]['scan start time'] > psm.RT_exp:
                    j -= 1
                basemz = spectra[j]
                I = []
                Ip = 0
                for idx, mz in enumerate(spectra[j]['m/z array']):
                    if abs(mz - (
                        float(psm.mass_exp + 1.007825 * psm.pcharge
                            ) / psm.pcharge)) <= 2:
                        if any(abs(mz - (float(
                            psm.mass_exp + k + 1.007825 * psm.pcharge
                            ) / psm.pcharge)) <= mz * 10e-5 for k in range(-2, 3)):
                            Ip += float(spectra[j]['intensity array'][idx])
                        else:
                            pass
                            #print abs(mz - (float(psm['precursor_neutral_mass']
                            #+ 1.007825 * psm['assumed_charge']) / psm['assumed_charge'])
                            #), psm['retention_time_sec'] / 60, spectra[j][
                            #'scanList']['scan'][0]['scan start time']
                        I.append(float(spectra[j]['intensity array'][idx]))
                PIF = Ip / sum(I) * 100
                psm.I = Ip
                psm.PIF = PIF

    spectra_dict = {}
    if len(argv) > 4 and argv[4].split('.')[-1] == 'mgf':
        fname = argv[4]
        spectra = mgf.read(fname)
        for spectrum in spectra:
#            print spectrum['params']['title'].split(',')[0]
            spectra_dict[spectrum['params']['title'].split(',')[0]] = spectrum['m/z array']

    print 'total number of PSMs = %d' % (len(peptides.peptideslist),)

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
            if protein.dbname not in valid_proteins and peptide.note == 'decoy': #### OR REPLACED AND
                protein.note = 'W'
                if peptide.note2 != 'tr' or peptide.note == 'decoy':
                    peptide.note2 = 'wr'
            else:
                true_prots.add(protein.dbname)
                protein.note = 'Valid'
                peptide.note2 = 'tr'

    for peptide in peptides.peptideslist:
        """
        try:
            spN = peptide.spectrum.split('Cmpd ')[1].split(',')[0]
        except:
            spN = str(int(peptide.spectrum.split('scan=')[1].split('"')[0]) - 1)
#            print spN
        """
        peptide.sumI = Fragment_intensities[peptide.start_scan]
        try:
            peptide.spectrum_mz = spectra_dict[peptide.spectrum.split(',')[0]]
        except:
            peptide.spectrum_mz = spectra_dict[
                    'Cmpd ' + peptide.spectrum.split(',')[0].split('.')[1]]


    for peptide in peptides.peptideslist:
        peptide.peptscore2 = pepts_dict[peptide.sequence]
        for protein in peptide.parentproteins:
            try:
                if peptide.protscore2 < float(
                        prots_dict[protein.dbname]
                        ) / protsL[protein.dbname] * 500:
                    peptide.protscore2 = float(
                            prots_dict[protein.dbname]
                            ) / protsL[protein.dbname] * 500
            except:
                print ('protein %s is missing in fasta, '
                        '5000 length is used for normalization' % protein.dbname)
                if peptide.protscore2 < float(prots_dict[protein.dbname]) / 5000 * 500:
                    peptide.protscore2 = float(prots_dict[protein.dbname]) / 5000 * 500
    copy_peptides = PeptideList()
    copy_peptides.peptideslist = list(peptides.peptideslist)
    copy_peptides.pepxml_type = peptides.pepxml_type
    copy_peptides.filter_evalue_new(FDR=FDR, useMP=False)

    print 'Default filtering:'
    PSMs_info(copy_peptides)

    descriptors = []
    descriptors.append(Descriptor(name='RT difference, min',
        formula='(peptide.RT_exp - peptide.RT_predicted)',
        group='A', binsize='auto'))
    descriptors.append(Descriptor(name='precursor mass difference, ppm',
        formula='peptide.mass_diff()', group='A', binsize='auto'))
    descriptors.append(Descriptor(name='missed cleavages',
        formula='peptide.num_missed_cleavages', group='A', binsize='1'))
#    descriptors.append(Descriptor(name='charge states',
#       formula='peptide.pcharge', group='A', binsize='1'))
    descriptors.append(Descriptor(name='potential modifications',
        formula="(peptide.modified_code.count('[') - peptide.modified_code.count('[160]'))",
        group='A', binsize='1'))
    descriptors.append(Descriptor(name='isotopes mass difference, Da',
        formula='round(peptide.massdiff, 0)', group='A', binsize='1'))
    descriptors.append(Descriptor(name='PSMs per protein',
        formula="peptide.protscore2", group='B'))
    descriptors.append(Descriptor(name='PSM count',
        formula='peptide.peptscore2', group='B'))
    descriptors.append(Descriptor(name='fragment mass tolerance, Da',
        formula='peptide.get_fragment_mt()', group='A', binsize = 'auto'))
#    descriptors.append(Descriptor(name='PIF', formula='peptide.PIF', group='B'))
#    descriptors.append(Descriptor(name='I', formula='peptide.I', group='B'))

    if 'RT difference, min' in [d.name for d in descriptors]:
        if RT_type == 'biolccc':
            copy_peptides.get_RC()
            print copy_peptides.RC
        #        copy_peptides.RC = pickle.load(open('RC.pickle'))
        #        pickle.dump(copy_peptides.RC, open('/home/mark/RC.pickle', 'w'))
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


    descriptors = prepare_hist(descriptors, copy_peptides)
    fig = plot_histograms(descriptors, peptides)

    if len(copy_peptides.peptideslist) > 50:
        jk = {}
        for peptide in peptides.peptideslist:
            jk = calc_peptscore(peptide, descriptors, jk)
        print jk
        j = len(peptides.peptideslist) - 1
        while j >= 0:
            if peptides.peptideslist[j].peptscore == 0:
                peptides.peptideslist.pop(j)
            j -= 1


    copy_peptides = PeptideList()
    copy_peptides.peptideslist = list(peptides.peptideslist)
    copy_peptides.pepxml_type = peptides.pepxml_type
    copy_peptides.filter_evalue_new(FDR=FDR, useMP=False)
    descriptors = prepare_hist(descriptors, copy_peptides)
    jk = {}
    for peptide in peptides.peptideslist:
        jk = calc_peptscore(peptide, descriptors, jk)
    print jk
#    plot_histograms(descriptors, copy_peptides)
    plot_MP(descriptors, peptides, fig)


main(inputfile, preFDR, FDR, RT_type, min_charge, max_charge)
