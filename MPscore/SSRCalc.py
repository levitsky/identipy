# -*- coding: utf-8 -*-
import os
import shelve
import mechanize
import urllib2
import re
import commands

DATABASE_FILENAME = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 'ssrcalc.db')


def remove_html_tags(data):
    p = re.compile(r'<.*?>')
    return p.sub('', data)


def calculate_RH(seq_list, pore_size=100, ion_pairing_agent='TFA', pH=2,
                proxy=''):
    # Find cached RHs in database.
    output = dict([(seq, None) for seq in seq_list])
    database = shelve.open(DATABASE_FILENAME, writeback=True)
    for seq in seq_list:
        if (seq in database
            and (pore_size, ion_pairing_agent, pH) in database[seq]):

            output[seq] = database[seq][(pore_size, ion_pairing_agent, pH)]
    remaining_seq = []
    for seq, RH in output.items():
        if RH is None:
            remaining_seq.append(seq)

    if not remaining_seq:
        return output

    # If there undefined RHs, obtain them at the SSRCalc site:
    if proxy:
        proxy_handler = urllib2.ProxyHandler({'http': proxy})
        opener = urllib2.build_opener(proxy_handler)
        opener.addheaders = [
            ('Host', '2ip.ru\n'),
            ('User-Agent',
                'Mozilla/5.0 (Windows; U; Windows NT 5.1; ru; rv:1.8.0.2) '
                'Gecko/20060308 Firefox/1.5.0.2\n'),
            ('Accept',
                'text/xml,application/xml,application/xhtml+xml,text/html;'
                'q=0.9,text/plain;q=0.8,image/png,*/*;q=0.5\n'),
            ('Accept-Language', 'ru-ru,ru;q=0.8,en-us;q=0.5,en;q=0.3\n'),
            ('Accept-Charset', 'windows-1251,utf-8;q=0.7,*;q=0.7\n'),
            ('X-Forwarded-For', '44.55.66.77\n'),
            ('Pragma', 'no-cache\n'),
            ('Referer', 'http://www.test.com\n'),
            ('Keep-Alive', '500\n'),
            ('Connection', 'close\n'),
            ('Content-Type', 'application/x-www-form-urlencoded\r\n\r\n')]
        urllib2.install_opener(opener)
    request = urllib2.Request('http://hs2.proteome.ca/SSRCalc/SSRCalcX.html')
    response = urllib2.urlopen(request)
    forms = mechanize.ParseResponse(response, backwards_compat=False)
    response.close()
    form = forms[0]
    if ion_pairing_agent == 'FA':
        form['sver'] = ['ssrFA']
    elif pH == 10:
        form['sver'] = ['ssrXT']
    elif pore_size == 100:
        form['sver'] = ['ssr100']
    elif pore_size == 300:
        form['sver'] = ['ssr300']
    form['seqs'] = "\n".join(remaining_seq)
    result = urllib2.urlopen(form.click())
    result = result.read()

    processed_seq_re = re.compile(
        r'(?<=\<tr class\=\"bodyText\"\>\<td\>)\S+\n?\S+')
    processed_RH_re = re.compile(r'(\(\d+\)\<\/td\>\n\<td\>\s+)(-?\d+.?\d+)')

    processed_seq = processed_seq_re.findall(result)
    processed_RH = [float(rh[1]) for rh in processed_RH_re.findall(result)]
    processed_data = dict(zip(processed_seq, processed_RH))

    # Caching obtained data.
    for seq, RH in processed_data.items():
        seq = remove_html_tags(seq)
        entry = database.get(seq, {})
        entry[(pore_size, ion_pairing_agent, pH)] = RH
        database[seq] = entry
        output[seq] = RH

    database.close()

    return output


def set_RH_psm_list(psm_list, pore_size, ion_pairing_agent, pH, proxy='',
                    local_version=None):
    seq_list = [x['peptide'] for x in psm_list]
    if local_version:
        RHs = calculate_RH_local(seq_list, local_version)
    else:
        RHs = calculate_RH(seq_list, pore_size, ion_pairing_agent, pH, proxy)
    for i in range(len(psm_list)):
        psm_list[i]['RH_SSRCalc'] = RHs[i]
    return psm_list


def calculateRH_local(seq_list, version=None):
    if version == 'nocluster':
        output = commands.getoutput('perl ./SSRCalc3_nocluster.pl --alg 3.0'
            ' --seq %s --output tsv' % ('/'.join(seq_list)))
    elif version == 'noelectric':
        output = commands.getoutput('perl ./SSRCalc3_noelectric.pl --alg 3.0'
            ' --seq %s --output tsv' % ('/'.join(seq_list)))
    elif version == 'nohelix':
        output = commands.getoutput('perl ./SSRCalc3_nohelix.pl --alg 3.0'
            ' --seq %s --output tsv' % ('/'.join(seq_list)))
    elif version == 'noall':
        output = commands.getoutput('perl ./SSRCalc3_noall.pl --alg 3.0'
            ' --seq %s --output tsv' % ('/'.join(seq_list)))
    else:
        output = commands.getoutput('perl ./SSRCalc3.pl --alg 3.0'
            ' --seq %s --output tsv' % ('/'.join(seq_list)))
    psm_list = []
    for string in output.split('\n'):
        splitted_string = string.split()
        psm_list.append((splitted_string[0], float(splitted_string[2])))
    return psm_list

