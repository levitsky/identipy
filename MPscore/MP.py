from MPlib import PeptideList, Descriptor, get_settings
from sys import argv
from math import log
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import percentileofscore, scoreatpercentile
from os import path
import pickle
from pyteomics import mzml, fasta, auxiliary, mgf, parser
import operator

inputfile = str(argv[1])

protsC = {}
protsL = {}


def find_optimal_xy(descriptors):
    x, y = 1, 1
    while x * y < len(descriptors) + 1:
        if x > y:
            y += 1
        else:
            x += 1
    return x, y

def PSMs_info(peptides, valid_proteins, prots_dict=None, printresults=True, tofile=False):
    proteins_groups = []
    full_sequences = set()
    tostay = set()

    for peptide in peptides.peptideslist:
        full_sequences.add(peptide.sequence)
    while full_sequences:
        prots = dict()
        prots_peptides = dict()
        peptides_added = set()
        true_prots = set()
    #    unique = set()
        todel = set()
        Total_prots = set()
        Total_prots_2 = set()
        for peptide in peptides.peptideslist:
            if peptide.sequence in full_sequences:
                for protein in peptide.parentproteins:
                    try:
                        prots[protein.dbname]['PSMs'] += 1
                        prots[protein.dbname]['sumI'] += peptide.sumI
                    except:
                        prots[protein.dbname] = dict()
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
        for protein in sorted(prots.items(), key=lambda x: x[1]['Peptides'], reverse=True):
            for peptide in peptides.peptideslist:
                if peptide.sequence in full_sequences:
                    if protein[0] in [p.dbname for p in peptide.parentproteins]:
                        full_sequences.remove(peptide.sequence)
            tostay.add(protein[0])
            break
#        print len(full_sequences), len(tostay)
#    for group in proteins_groups:
#        if any([prots[dbname]['Peptides'] > 1 for dbname in group]):
#            print group

    prots = dict()
    prots_peptides = dict()
    peptides_added = set()
    true_prots = set()
#    unique = set()
    todel = set()
    Total_prots = set()
    Total_prots_2 = set()
    for peptide in peptides.peptideslist:
        peptide.sequence = peptide.sequence.replace('|','')
        flag = 1
        for group in proteins_groups:
            if any(protein.dbname in group for protein in peptide.parentproteins):
                for protein in peptide.parentproteins:
                    group.add(protein.dbname)
                flag = 0
                break
        if flag:
            proteins_groups.append(set(protein.dbname for protein in peptide.parentproteins))

        for protein in peptide.parentproteins:
            if protein.dbname in Total_prots:
                Total_prots_2.add(protein.dbname)
            else:
                Total_prots.add(protein.dbname)
            try:
                prots[protein.dbname]['PSMs'] += 1
                prots[protein.dbname]['sumI'] += peptide.sumI
            except:
                prots[protein.dbname] = dict()
                prots[protein.dbname]['PSMs'] = 1
                prots[protein.dbname]['sumI'] = peptide.sumI
            if protein.dbname in valid_proteins and peptide.note != 'decoy':
                true_prots.add(protein.dbname)
        if peptide.sequence not in peptides_added:
            for protein in peptide.parentproteins:
                try:
                    prots[protein.dbname]['Peptides'] += 1
                except:
                    prots[protein.dbname]['Peptides'] = 1
        peptides_added.add(peptide.sequence)

    for dbname in list(prots.keys()):
        if dbname not in tostay:
            del prots[dbname]
    if tofile:
        output_proteins = open('%s/%s_proteins.csv' % (path.dirname(path.realpath(inputfile)), path.splitext(path.splitext(path.basename(inputfile))[0])[0]), 'w')
        output_peptides_detailed = open('%s/%s_peptides.csv' % (path.dirname(path.realpath(inputfile)), path.splitext(path.splitext(path.basename(inputfile))[0])[0]), 'w')
        output_peptides_detailed.write('sequence\tmodified_sequence\te-value\tMPscore\tspectrum_title\tproteins\n')
        if protsC:
            output_proteins_valid = open('%s/Results_new_%s_proteins_valid.csv' % (path.dirname(path.realpath(inputfile)), path.splitext(path.splitext(path.basename(inputfile))[0])[0]), 'w')
            temp_data = []
        for k, v in prots.items():
            if protsC and k in valid_proteins:
                output_proteins_valid.write('%s,%s,%s,%s,%s\n' % (k, v['PSMs'], v['Peptides'], v['sumI'], protsC[k]))
                temp_data.append([float(v['sumI']), protsC[k]])
            if int(v['Peptides']) > 1:
                output_proteins.write('%s\t%s\t%s\t%s\t%s\n' % (k, v['PSMs'], v['Peptides'], v['sumI'], ('+' if k in valid_proteins else '-')))

        peptides_best = dict()
        peptides_best_sp = dict()
        for peptide in peptides.peptideslist:
            if peptide.sequence in peptides_best.keys():
                if peptide.evalue < peptides_best[peptide.sequence]:
                    peptides_best[peptide.sequence] = peptide.evalue
                    peptides_best_sp[peptide.sequence] = peptide.spectrum
            else:
                peptides_best[peptide.sequence] = peptide.evalue
                peptides_best_sp[peptide.sequence] = peptide.spectrum
        for peptide in peptides.peptideslist:
            if peptide.spectrum == peptides_best_sp[peptide.sequence]:
                output_peptides_detailed.write('%s\t%s\t%s\t%s\t%s\t' % (peptide.sequence, peptide.modified_code, peptide.evalue, peptide.peptscore, peptide.spectrum))
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
            print 'conc: ', auxiliary.linear_regression([x[0] for x in temp_data], [x[1] for x in temp_data])

        output_proteins.close()
    if printresults:
        print 'True PSMs: %s' % (len([1 for x in peptides.peptideslist if x.note2 == 'tr']), )
        print 'True2 PSMs: %s' % (len([1 for x in peptides.peptideslist if x.note3]), )
        print 'Peptides: %d' % (len(set(p.sequence for p in peptides.peptideslist)))
#        print 'Wrong Peptides= %s' % (len(set([x.sequence for x in peptides.peptideslist if not x.note3])), )
        print 'Protein groups: %s' % (len(prots.values()))
        print 'Protein groups with >= 2 peptides: %s' % (len([v for v in prots.values() if v['Peptides'] >= 2]))
    #    print 'Total_prots >=2 = %s' % (len(Total_prots_2), )
        if valid_proteins:
            print 'True Prots = %s' % (len(true_prots))
            print 'Real FDR = %s' % (100 * float(len([1 for x in peptides.peptideslist if not x.note3])) / len(peptides.peptideslist) )
        print '\n'
    return (len([1 for x in peptides.peptideslist if x.note2 == 'tr']), len(set(p.sequence for p in peptides.peptideslist)), len([v for v in prots.values() if v['Peptides'] > 1]))

def plot_histograms(descriptors, peptides, FDR):
    fig = plt.figure(figsize=(16, 12))
    ox, oy = find_optimal_xy(descriptors)
    copy_peptides = PeptideList()
    copy_peptides.peptideslist = list(peptides.peptideslist)
    copy_peptides.pepxml_type = peptides.pepxml_type
    copy_peptides.filter_evalue_new(FDR1=FDR, useMP=False)
    for idx, descriptor in enumerate(descriptors):
        ax = fig.add_subplot(ox, oy, idx + 1)
        array_wrong = [eval(descriptor.formula) for peptide in peptides.peptideslist if peptide.note2 == 'wr']# and peptide.note == 'decoy']
        array_valid = [eval(descriptor.formula) for peptide in peptides.peptideslist if peptide.note2 == 'tr']
        if descriptor.group == 'B':
            array_wrong = np.log10(array_wrong)
            array_valid = np.log10(array_valid)
        binsize = float(descriptor.get_binsize(copy_peptides))
        if binsize < float(max(np.append(array_wrong, array_valid)) - min(np.append(array_wrong, array_valid))) / 400:
            binsize = float(max(np.append(array_wrong, array_valid)) - min(np.append(array_wrong, array_valid))) / 400
        lbin_s = scoreatpercentile(np.append(array_wrong, array_valid), 1.0)
        lbin = min(np.append(array_wrong, array_valid))
        if lbin_s and abs((lbin - lbin_s) / lbin_s) > 1.0:
            lbin = lbin_s * 1.05
        rbin_s = scoreatpercentile(np.append(array_wrong, array_valid), 99.0)
        rbin = max(np.append(array_wrong, array_valid))
        if rbin_s and abs((rbin - rbin_s) / rbin_s) > 1.0:
            rbin = rbin_s * 1.05
        rbin += 1.5 * binsize
        if descriptor.name == 'fragment mass tolerance, Da':
            rbin = 0.05
        H1, _ = np.histogram(array_wrong, bins=np.arange(lbin, rbin, binsize))
        H2, _ = np.histogram(array_valid, bins=np.arange(lbin, rbin, binsize))
        if descriptor.group == 'B':
            H3, _ = np.histogram([np.log10(eval(descriptor.formula)) for peptide in copy_peptides.peptideslist], bins=np.arange(lbin, rbin, binsize))
        else:
            H3, _ = np.histogram([eval(descriptor.formula) for peptide in copy_peptides.peptideslist], bins=np.arange(lbin, rbin, binsize))
        if descriptor.name in ['precursor mass difference, ppm']:
            print 'MEDIAN precursor mass difference of top PSMs=%s ppm' % (np.median([eval(descriptor.formula) for peptide in copy_peptides.peptideslist]), )
        if descriptor.group == 'B':
            H1 = H1.clip(1)
            H2 = H2.clip(1)
            H3 = H3.clip(1)
            H1 = np.log10(H1)
            H2 = np.log10(H2)
            H3 = np.log10(H3)
        ind = np.arange(lbin, rbin - (1 + 1e-15) * binsize, binsize)
        width = binsize
        plt.bar(ind, H1, width, color='red', alpha=0.5)
        plt.bar(ind, H2, width, color='green', alpha=0.5)
        plt.bar(ind, H3, width, color='black', alpha=0.5)

        if descriptor.name in ['missed cleavages', 'charge states', 'potential modifications', 'isotopes mass difference, Da']:
            plt.xticks(range(1, 6))
        ax.set_xlabel(descriptor.name)
    return fig


def plot_MP(descriptors, peptides, fig, FDR, FDR_new, prots_dict, valid_proteins, threshold0=False):
    ox, oy = find_optimal_xy(descriptors)
    copy_peptides = PeptideList()
    copy_peptides.peptideslist = list(peptides.peptideslist)
    copy_peptides.pepxml_type = peptides.pepxml_type
    threshold1, threshold2 = copy_peptides.filter_evalue_new(FDR1=FDR, FDR2=FDR_new, useMP=True)
    threshold1 = -log(threshold1)
    try:
        threshold2 = log(threshold2)
    except:
        pass

    zero_peptscore = log(min([p.peptscore for p in peptides.peptideslist if p.peptscore != 0])) - 1
    zero_evalue = -log(min([p.evalue for p in peptides.peptideslist if p.evalue != 0])) + 1
    PSMs_wrong = [[(-log(pept.evalue) if pept.evalue != 0 else zero_evalue), (log(pept.peptscore) if pept.peptscore != 0 else zero_peptscore)] for pept in peptides.peptideslist if pept.note2 == 'wr']
    PSMs_true = [[(-log(pept.evalue) if pept.evalue != 0 else zero_evalue), (log(pept.peptscore) if pept.peptscore != 0 else zero_peptscore)] for pept in peptides.peptideslist if pept.note2 == 'tr']

    print 'MP filtering:'
    PSMs_info(copy_peptides, valid_proteins, prots_dict=prots_dict, tofile=True)
    print 'Without filtering, after removing outliers:'
    PSMs_info(peptides, valid_proteins)

    ax = fig.add_subplot(ox, oy, ox * oy)
    ax.plot([x[0] for x in PSMs_wrong], [x[1] for x in PSMs_wrong], 'o', markersize=2, color='red')
    ax.plot([x[0] for x in PSMs_true], [x[1] for x in PSMs_true], 'o', markersize=2, color='blue')
    ax.axvline(threshold1, color='green')
    if threshold2:
        ax.axhline(threshold2, color='green')
    if threshold0:
        threshold0 = -log(threshold0)
        ax.axvline(threshold0, color='red')
    ax.set_ylim(min(x[1] for x in PSMs_true) - 1, max(x[1] for x in PSMs_true) + 1)
    fig.tight_layout()
    plt.savefig('%s/Results_new_%s.png' % (path.dirname(path.realpath(inputfile)), path.splitext(path.splitext(path.basename(inputfile))[0])[0]))

def prepare_hist(descriptors, copy_peptides, first=False):
    for descriptor in descriptors:
        descriptor.array = descriptor.get_array(copy_peptides)
        if descriptor.group == 'A':
            binsize = float(descriptor.get_binsize(copy_peptides))
            if binsize < float(max(descriptor.get_array(copy_peptides)) - min(descriptor.get_array(copy_peptides))) / 400:
                    binsize = float(max(descriptor.get_array(copy_peptides)) - min(descriptor.get_array(copy_peptides))) / 400

            lbin, rbin = min(descriptor.array), max(descriptor.array) + 1.5 * binsize
            if lbin == rbin:
                rbin = lbin + binsize
            descriptor.hist = np.histogram(descriptor.array, bins=np.arange(lbin, rbin + binsize, binsize))

            if descriptor.name.startswith('potential modifications'):
                descriptor.hist[0][0] = max(descriptor.hist[0][1:])
                if not descriptor.hist[0][0]:
                    descriptor.hist[0][0] = 1
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
                            peptide.peptscore *= float(descriptor.hist[0][j]) / sum(descriptor.hist[0])
                            break
                        j -= 1
        elif descriptor.group == 'B':
            peptide.peptscore *= percentileofscore(descriptor.array, eval(descriptor.formula)) / 100

        if flag and not peptide.peptscore:
            try:
                jk[descriptor.name] += 1
            except:
                jk[descriptor.name] = 1

    return jk



def main(inputfile):

    print 'inputfile = %s' % (inputfile, )
    basename = path.splitext(path.splitext(argv[1])[0])[0]
    txmlfile = None
    mzmlfile = None
    mgffile = None
    fastafile = None
    configfile = None

    for arg in argv:
        if path.splitext(arg)[-1] == '.xml':
            if path.splitext(path.splitext(arg)[0])[-1] == '.t':
                txmlfile = arg
        elif path.splitext(arg)[-1] == '.mzml':
            mzmlfile = arg
        elif path.splitext(arg)[-1] == '.mgf':
            mgffile = arg
        elif path.splitext(arg)[-1] == '.fasta':
            fastafile = arg
        elif path.splitext(arg)[-1] == '.cfg':
            configfile = arg

    if not txmlfile:
        if path.isfile(basename + path.extsep + 't' + path.extsep + 'xml'):
            txmlfile = basename + path.extsep + 't' + path.extsep + 'xml'
    if not mzmlfile:
        if path.isfile(basename + path.extsep + 'mzml'):
            mzmlfile = basename + path.extsep + 'mzml'
    if not mgffile:
        if path.isfile(basename + path.extsep + 'mgf'):
            mgffile = basename + path.extsep + 'mgf'

    if configfile:
        settings = get_settings(configfile)
    else:
        settings = get_settings('default.cfg')

    valid_proteins = []
    valid_proteins_input = settings.get('options', 'valid proteins')
    if path.isfile(valid_proteins_input):
        for line in open(valid_proteins_input, 'r'):
            if len(line.strip().split()) > 1:
                dbname, conc = line.strip().split()
                valid_proteins.append(dbname)
                protsC[dbname] = float(conc)
            else:
                dbname = line.strip().split()[0]
                valid_proteins.append(dbname)
    else:
        valid_proteins = []



    line = inputfile
    search_space_peptides, decoy_PSMs_all = set(), 0.
    numprots, numpeptides, numprots_true, numpeptides_true = 0, 0, 0, 0
    FDR = settings.getfloat('options', 'FDR')

    proteases = [x.strip() for x in settings.get('missed cleavages', 'protease1').split(',')]
    proteases.extend([x.strip() for x in settings.get('missed cleavages', 'protease2').split(',') if x.strip()])
    expasy = '|'.join((parser.expasy_rules[protease] if protease in parser.expasy_rules.keys() else protease for protease in proteases))

    RT_type = settings.get('retention time', 'model')

    min_charge = settings.getint('charges', 'min charge')
    max_charge = settings.getint('charges', 'max charge')

    peptides = PeptideList(settings)
    peptides.get_from_pepxmlfile(line, min_charge=min_charge, max_charge=max_charge, max_rank=1)
    missed_modifications = set()
    for peptide in peptides.peptideslist:
        for mod in peptide.missed_modifications:
            missed_modifications.add(mod)
    for mod in missed_modifications:
        print 'modification with mass = %s is missed in parameters' % (mod)

    Fragment_intensities = {}
    spectra_dict = dict()

    if fastafile:
        for x in fasta.read(fastafile):
            if not any(x[0].startswith(tag) for tag in ['sp', 'tr', 'DECOY_sp', 'DECOY_tr']):
                if any(tag in x[0] for tag in ['SWISS-PROT:', 'TREMBL:']):
                    dbname = x[0].split(' ')[0]
                else:
                    dbname = x[0]
                protsL[dbname] = len(x[1])
            else:
                protsL[x[0].split('|')[1]] = len(x[1])

            numprots += 1
#            for peptide in parser.cleave(x[1], expasy, 2):
#                numpeptides += 1
            search_space_peptides.update(parser.cleave(x[1], expasy, 2))
        numpeptides = len(search_space_peptides)

    if txmlfile:
        for x in open(txmlfile, 'r').readlines():
            if x.startswith('<group id='):
                Fragment_intensities[int(x.split('<group id="')[1].split('"')[0])] = 10**float(x.split('sumI="')[1].split('"')[0])

    if mzmlfile:
        isolation_window = settings.getfloat('precursor ion fraction', 'isolation window')
        mass_acc = settings.getfloat('precursor ion fraction', 'mass accuracy')
        spectra = [x for x in mzml.read(mzmlfile) if x['ms level'] == 1]
        for psm in peptides.peptideslist:
            j = len(spectra) - 1
            while spectra[j]['scanList']['scan'][0]['scan start time'] > psm.RT_exp:
                j -= 1
            basemz = spectra[j]
            I = []
            Ip = 0
            for idx, mz in enumerate(spectra[j]['m/z array']):
                if abs(mz - (float(psm.mass_exp + 1.007825 * psm.pcharge) / psm.pcharge)) <= isolation_window:
                    if any(abs(mz - (float(psm.mass_exp + k + 1.007825 * psm.pcharge) / psm.pcharge)) <= mz * mass_acc * 1e-6 for k in [-2, -1, 0, 1, 2]):
                        Ip += float(spectra[j]['intensity array'][idx])
                    I.append(float(spectra[j]['intensity array'][idx]))
            PIF = Ip / sum(I) * 100
            psm.I = Ip
            psm.PIF = PIF

    if mgffile:
        spectra = mgf.read(mgffile)
        spectrum = spectra.next()
        if len(spectrum['params']['title'].split(',')) > 1:
            spectra_name_type = 'Bruker mgf'
            spectra_dict[int(spectrum['params']['title'].split(',')[0].split('Cmpd ')[1])] = spectrum['m/z array']
        elif len(spectrum['params']['title'].split('.')) > 1 and spectrum['params']['title'].split('.')[1] == spectrum['params']['title'].split('.')[2]:
            spectra_name_type = 'TPP'
            spectra_dict[int(spectrum['params']['title'].split('.')[1])] = spectrum['m/z array']
        else:
            print 'Unknown spectra name type in mgf file, you have to turn off fragment mass tolerance error calculation!'
        for spectrum in spectra:
            if spectra_name_type == 'Bruker mgf':
                spectra_dict[int(spectrum['params']['title'].split(',')[0].split('Cmpd ')[1])] = spectrum['m/z array']
            elif spectra_name_type == 'TPP':
                spectra_dict[int(spectrum['params']['title'].split('.')[1])] = spectrum['m/z array']

    print 'total number of PSMs = %d' % (len(peptides.peptideslist),)

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
                if protein.dbname in valid_proteins:
                    peptide.note3 = 'valid'
                true_prots.add(protein.dbname)
                protein.note = 'Valid'
                peptide.note2 = 'tr'

    for peptide in peptides.peptideslist:
        if peptide.note == 'decoy':
            decoy_PSMs_all += 1
        try:
            peptide.sumI = Fragment_intensities[peptide.start_scan]
        except:
            peptide.sumI = 0

    if settings.getint('descriptors', 'fragment mass tolerance, Da'):
        for peptide in peptides.peptideslist:
            """
            try:
                spN = peptide.spectrum.split('Cmpd ')[1].split(',')[0]
            except:
                spN = str(int(peptide.spectrum.split('scan=')[1].split('"')[0]) - 1)
    #            print spN
            """
            if spectra_dict:
    #            try:
                if spectra_name_type == 'TPP':
                    try:
                        peptide.spectrum_mz = spectra_dict[int(peptide.spectrum.split('.')[1])]
                    except:
                        spectra_name_type = 'Bruker mgf'
                        peptide.spectrum_mz = spectra_dict[int(peptide.spectrum.split(',')[0].split('Cmpd ')[1])]
    #            except: 
                elif spectra_name_type == 'Bruker mgf':
    #                print peptide.spectrum
                    try:   
                        peptide.spectrum_mz = spectra_dict[int(peptide.spectrum.split(',')[0].split('Cmpd ')[1])]
                    except:
                        spectra_name_type = 'TPP'
                        peptide.spectrum_mz = spectra_dict[int(peptide.spectrum.split('.')[1])]
            else:
                print 'mgf file is missing'
            peptide.get_fragment_mt(peptides.settings)
            

    for peptide in peptides.peptideslist:
        peptide.peptscore2 = pepts_dict[peptide.sequence]
        for protein in peptide.parentproteins:
            try:
                if peptide.protscore2 < float(prots_dict[protein.dbname]) / protsL[protein.dbname] * 500:
                    peptide.protscore2 = float(prots_dict[protein.dbname]) / protsL[protein.dbname] * 500
            except:
                print 'protein %s is missed in fasta, 5000 length is used for normalization' % (protein.dbname, )
                if peptide.protscore2 < float(prots_dict[protein.dbname]) / 5000 * 500:
                    peptide.protscore2 = float(prots_dict[protein.dbname]) / 5000 * 500
    copy_peptides = PeptideList()
    copy_peptides.peptideslist = list(peptides.peptideslist)
    copy_peptides.pepxml_type = peptides.pepxml_type
    threshold0, _ = copy_peptides.filter_evalue_new(FDR1=FDR, useMP=False)

    print 'Default filtering:'
    PSMs_info(copy_peptides, valid_proteins)

    descriptors = []
    dname = 'RT difference, min'
    if settings.getint('descriptors', dname):
        descriptors.append(Descriptor(name=dname, formula='(peptide.RT_exp - peptide.RT_predicted)', group='A', binsize='auto'))
    dname = 'precursor mass difference, ppm'
    if settings.getint('descriptors', dname):
        descriptors.append(Descriptor(name=dname, formula='peptide.mass_diff()', group='A', binsize='auto'))
    dname = 'missed cleavages, protease 1'
    if settings.getint('descriptors', dname.split(',')[0]):
        protease1 = [x.strip() for x in settings.get('missed cleavages', 'protease1').split(',')]
        expasy1 = '|'.join((parser.expasy_rules[protease] if protease in parser.expasy_rules.keys() else protease for protease in protease1))
#        descriptors.append(Descriptor(name=dname, formula='peptide.get_missed_cleavages(%s)' % ('"'+protease+'"', ), group='A', binsize=1))
        descriptors.append(Descriptor(name=dname, formula='peptide.get_missed_cleavages(%s)' % ('"'+expasy1+'"', ), group='A', binsize=1))
    dname = 'missed cleavages, protease 2'
    if settings.getint('descriptors', dname.split(',')[0]):
        protease2 = [x.strip() for x in settings.get('missed cleavages', 'protease2').split(',')]
        expasy2 = '|'.join((parser.expasy_rules[protease] if protease in parser.expasy_rules.keys() else protease for protease in protease2))
        if protease2[0]:
            print protease2
            print expasy2
#        descriptors.append(Descriptor(name=dname, formula='peptide.get_missed_cleavages(%s)' % ('"'+protease+'"', ), group='A', binsize=1))
            descriptors.append(Descriptor(name=dname, formula='peptide.get_missed_cleavages(%s)' % ('"'+expasy2+'"', ), group='A', binsize=1))

    dname = 'charge states'
    if settings.getint('descriptors', dname):
        descriptors.append(Descriptor(name=dname, formula='peptide.pcharge', group='A', binsize='1'))
    dname = 'potential modifications'
    if settings.getint('descriptors', dname):
        labeldict = dict()
        temp = settings.get('modifications', 'variable')
        if temp:
            for mod in temp.split(', '):
                descriptors.append(Descriptor(name='%s, %s' % (dname, mod), formula="peptide.count_modifications_new(%s)" % ('"'+mod+'"', ), group='A', binsize='1'))
    dname = 'isotopes mass difference, Da'
    if settings.getint('descriptors', dname):
        descriptors.append(Descriptor(name=dname, formula='round(peptide.massdiff, 0)', group='A', binsize='1'))
    dname = 'PSMs per protein'
    if settings.getint('descriptors', dname):
        descriptors.append(Descriptor(name=dname, formula="peptide.protscore2", group='B'))
    dname = 'PSM count'
    if settings.getint('descriptors', dname):
        descriptors.append(Descriptor(name=dname, formula='peptide.peptscore2', group='B'))
    dname = 'fragment mass tolerance, Da'
    if settings.getint('descriptors', dname):
        if not spectra_dict:
            print 'mgf file is missing. Fragment mass tolerance descriptor is turned off.'
        else:
            descriptors.append(Descriptor(name=dname, formula='peptide.get_fragment_mt()', group='A', binsize = 'auto'))
    dname = 'PIF'
    if settings.getint('descriptors', dname):
        descriptors.append(Descriptor(name=dname, formula='peptide.PIF', group='B'))

    if 'RT difference, min' in [d.name for d in descriptors]:
        if RT_type == 'achrom':
            copy_peptides.get_RC()
            peptides.RC = copy_peptides.RC
            copy_peptides.calc_RT(RTtype=RT_type)
            print copy_peptides.get_calibrate_coeff()
            peptides.calc_RT(RTtype=RT_type)
        else:
            copy_peptides.filter_modifications()
            copy_peptides.calc_RT(RTtype=RT_type)
            calibrate_coeff = copy_peptides.get_calibrate_coeff()
            copy_peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)
            copy_peptides.filter_RT(RT_tolerance=3 * calibrate_coeff[3])
            copy_peptides.calc_RT(RTtype=RT_type)
            calibrate_coeff = copy_peptides.get_calibrate_coeff()
            print calibrate_coeff
            peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)
            copy_peptides.peptideslist = list(peptides.peptideslist)
            copy_peptides.filter_evalue_new(FDR1=FDR, useMP=False)
            copy_peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)


    descriptors = prepare_hist(descriptors, copy_peptides, first=False)
    fig = plot_histograms(descriptors, peptides, FDR)

    if len(copy_peptides.peptideslist) > 100:
        jk = dict()
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
    copy_peptides.filter_evalue_new(FDR1=FDR, useMP=False)
    descriptors = prepare_hist(descriptors, copy_peptides)
    jk = dict()
    for peptide in peptides.peptideslist:
        jk = calc_peptscore(peptide, descriptors, jk)
#    print jk
#    plot_histograms(descriptors, copy_peptides)
    k_temp = []
    while len(k_temp) < 3 or k_temp[-1] != k_temp[-3]:
        copy_peptides = PeptideList()
        copy_peptides.peptideslist = list(peptides.peptideslist)
        copy_peptides.pepxml_type = peptides.pepxml_type
        if not k_temp:
            copy_peptides.filter_evalue_new(FDR1=FDR, useMP=False)
        else:
            FDR_new = (((FDR / 100 - float(decoy_PSMs_all) / numPSMs_true * k)) / (1 - k) ) * 100
            copy_peptides.filter_evalue_new(FDR1=FDR, FDR2=FDR_new, useMP=True)

        numpeptides_true = set()
        numPSMs_true = 0.
        for peptide in copy_peptides.peptideslist:
            numpeptides_true.add(peptide.sequence)
            numPSMs_true += 1
        k = float(len(numpeptides_true)) / numpeptides
        FDR_new = (((FDR / 100 - float(decoy_PSMs_all) / numPSMs_true * k)) / (1 - k) ) * 100
        k_temp.append(float(k))
    print k, 'k factor'
    print FDR_new, 'FDR_new'
    plot_MP(descriptors, peptides, fig, FDR, FDR_new, prots_dict, valid_proteins, threshold0)

main(inputfile)
