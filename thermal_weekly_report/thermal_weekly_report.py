
import sys
import numpy as np
import jinja2 as ja
# from BeautifulSoup import BeautifulSoup as bs
from os.path import join as pathjoin
from os.path import expanduser
import cPickle as pickle
from shutil import copyfile
import argparse
import os

import Ska.engarchive.fetch_eng as fetch
from Chandra.Time import DateTime

home = expanduser("~")
sys.path.append(home + '/AXAFLIB/pylimmon/')
sys.path.append(home + '/AXAFLIB/pyeclipse/')

import pyeclipse as ecl
import pylimmon


# Patch BeautifulSoup prettify to increase indent
# orig_prettify = bs.prettify
# r = re.compile(r'^(\s*)', re.MULTILINE)
# def prettify(s, indent_width=4):
#     ret = r.sub(r'\1' * indent_width, orig_prettify(s))
#     return ret
# bs.prettify = prettify


def get_average_tel_power(t1, t2):
    """ Get average telescope power used during period

    :param t1: start of period
    :param t2: end of period
    :return: mean total (PTOTAL) telescope power for period
    """
    p = fetch.Msid('DP_PTOTAL', DateTime(t1).date, DateTime(t2).date,
                   stat='daily')
    return np.mean(p.vals)


def writeEclipseText(eclipse, ecllist):

    ecltext = []
    for eclnum in ecllist:
        if eclipse[eclnum].has_key('umbra'):
            tstart = eclipse[eclnum]['entrancepenumbra']['Start Time']
            tstop = eclipse[eclnum]['exitpenumbra']['Stop Time']
            text = "Umbral Eclipse: Start Time = %s, End Time = %s<br>\n"\
                   % (tstart, tstop)
        else:
            tstart = eclipse[eclnum]['entrancepenumbra']['Start Time']
            tstop = eclipse[eclnum]['entrancepenumbra']['Stop Time']
            text = "Penumbral Eclipse: Start Time = %s, End Time = %s<br>\n"\
                   % (tstart, tstop)

        ecltext.append(text)

    return ''.join(ecltext)


def get_eclipse_text(eclfile, t1, t2):

    #----------------------------------
    # Import eclipse data
    eclipse = ecl.read_eclipse_file(eclfile)
    # eclipse = ecl.read_eclipse_file('./ECLIPSE_HISTORY_2014.txt')

    #----------------------------------
    # Look for eclipses that fall in the current week
    ecllist = []
    for k in eclipse['eclipse_nums']:
        tstart = eclipse[k]['entrancepenumbra']['startsec']
        if eclipse[k].has_key('exitpenumbra'):
            # This is used for a normal eclipse
            tstop = eclipse[k]['exitpenumbra']['stopsec']
        else:
            # This is used for a penumbral eclipse
            tstop = eclipse[k]['entrancepenumbra']['stopsec']

        if (tstart >= DateTime(t1).secs) & (tstop <= DateTime(t2).secs):
            ecllist.append(k)

    if ecllist:
        eclipselines = writeEclipseText(eclipse, ecllist)
    else:
        eclipselines = 'None.'

    return eclipselines


def gen_thermal_checklist(filename='thermlist.csv'):

    def read_therm_checklist(filename):
        def splitline(line):
            return line.strip().split(',')

        with open(filename, 'r') as fid:
            lines = fid.readlines()

        header = splitline(lines.pop(0))
        thermlist = [splitline(line) for line in lines]

        return thermlist, header

    thermlist, header = read_therm_checklist(filename)

    thermdict = {}
    notinarchive = []
    missing = []
    for i, line in enumerate(thermlist):
        ok = False
        greta_msid = line[0].strip().lower()
        ska_msid = line[1].strip().lower()
        try:
            # Assume the local archive is not more than 30 days out of date
            data = fetch.Msid(ska_msid, DateTime().secs - 3600 * 24 * 30, stat='daily')
            ok = True
        except IOError as e:
            missing.append(thermlist[i])
        except ValueError as e:
            notinarchive.append(thermlist[i])

        if ok:
            thermdict[ska_msid] = {
                'greta_msid': greta_msid, 'owner': line[2].strip(), 'description':line[3].strip()}
            if data.state_codes:
                thermdict[ska_msid]['type'] = 'expst'
            else:
                thermdict[ska_msid]['type'] = 'limit'

    with open('thermalmsiddata.pkl', 'w') as fid:
        pickle.dump((thermdict, missing, notinarchive), fid, protocol=2)


def check_violations(thermdict, t1, t2):
    t1 = DateTime(t1).date
    t2 = DateTime(t2).date

    allviolations = {}
    missingmsids = []
    checkedmsids = []
    for key in thermdict.keys():
        greta_msid = thermdict[key]['greta_msid']
        try:
            if thermdict[key]['type'] == 'limit':
                if "wide" in greta_msid.lower():
                    violations = handle_widerange_cases(key, t1, t2, greta_msid)
                    checkedmsids.append(key)
                else:
                    violations = pylimmon.check_limit_msid(key, t1, t2, greta_msid=greta_msid)
                    checkedmsids.append(key)
            elif thermdict[key]['type'] == 'expst':
                violations = pylimmon.check_state_msid(key, t1, t2, greta_msid=greta_msid)
                checkedmsids.append(key)

            if violations['any']:
                allviolations[key] = violations
#                 allviolations['any'] = True

        except IndexError:
            print('{} not in DB'.format(key))
            missingmsids.append(key)

    return allviolations, missingmsids, checkedmsids


def handle_widerange_cases(key, t1, t2, greta_msid):
    if DateTime(t2).secs <= DateTime('2014:342:16:30:00').secs:
        violations = pylimmon.check_limit_msid(key, t1, t2, greta_msid=key)
    elif DateTime(t1).secs >= DateTime('2014:342:16:33:00').secs:
        violations = pylimmon.check_limit_msid(key, t1, t2, greta_msid=greta_msid)
    else:
        t2_a = np.min((DateTime(t2).secs, DateTime('2014:342:16:30:00').secs))
        violations = pylimmon.check_limit_msid(key, t1, t2_a, greta_msid=key)
        t1_b = np.min((DateTime(t2).secs, DateTime('2014:342:16:33:00').secs))
        violations_b = pylimmon.check_limit_msid(key, t1_b, t2, greta_msid=greta_msid)

        for limittype in ['warning_low', 'caution_low', 'caution_high', 'warning_high']:
            if limittype in violations_b.keys():
                if limittype in violations.keys():
                    violations[limittype]['extrema'].extend(violations_b[limittype]['extrema'])
                    violations[limittype]['times '].extend(violations_b[limittype]['times '])
                    violations[limittype]['timespans '].extend(
                        violations_b[limittype]['timespans '])
                else:
                    violations[limittype] = violations_b[limittype]
                    violations['any'] = True

    return violations


def add_violation_info(allviolations):
    for key in allviolations.keys():
        try:
            t1 = DateTime().secs
            data = fetch.Msid(key, t1 - 3600 * 24 * 30, t1, stat='daily')
            desc = data.tdb.technical_name
        except:
            desc = 'No Description in TDB'
        allviolations[key]['description'] = desc

        if 'expst' in allviolations[key]['type']:
            allviolations[key]['total_duration'] = sum(
                [b - a for a, b in allviolations[key]['timespans']]) / 3600.
            allviolations[key]['num_excursions'] = len(allviolations[key]['timespans'])
            allviolations[key]['expectation'] = '= {}'.format(allviolations[key]['limits'][0])
            allviolations[key]['observed'] = ', '.join(np.unique(allviolations[key]['extrema']))
            allviolations[key]['datespans'] = [(DateTime(a).date, DateTime(b).date)
                                               for a, b in allviolations[key]['timespans']]
        else:
            for limtype in ['warning_low', 'caution_low', 'caution_high', 'warning_high']:
                if limtype in allviolations[key].keys():
                    allviolations[key][limtype]['total_duration'] = sum([b - a for a, b in
                        allviolations[key][limtype]['timespans']]) / 3600.
                    if 'high' in limtype:
                        allviolations[key][limtype]['observed'] = np.max(
                            allviolations[key][limtype]['extrema'])
                        allviolations[key][limtype]['expectation'] = '< {}'.format(
                            allviolations[key][limtype]['limits'][0])
                    elif 'low' in limtype:
                        allviolations[key][limtype]['observed'] = np.min(
                            allviolations[key][limtype]['extrema'])
                        allviolations[key][limtype]['expectation'] = '> {}'.format(
                            allviolations[key][limtype]['limits'][0])

                    allviolations[key][limtype]['num_excursions'] = len(
                        allviolations[key][limtype]['timespans'])

                    allviolations[key][limtype]['datespans'] = [(DateTime(a).date, DateTime(b).date)
                        for a, b in allviolations[key][limtype]['timespans']]
    return allviolations


def check_limit_changes(t1, t2):

    db = pylimmon.open_sqlite_file()
    cursor = db.cursor()
    cursor.execute('''SELECT a.msid, a.setkey FROM limits AS a WHERE a.datesec>=? 
                      AND a.datesec <=? ''', [DateTime(t1).secs, DateTime(t2).secs])
    allchanges = cursor.fetchall()

    msid_sets = [(d[0], d[1]) for d in allchanges]
    msid_sets = set(msid_sets)

    changes = {}

    for msid, setval in list(msid_sets):

        try:
            if 'wide' in msid.lower():
                skamsid = msid[:-5]
            else:
                skamsid = msid

            data = fetch.Msid(skamsid, t1, DateTime(t1).secs + 3600, stat='5min')
            desc = data.tdb.technical_name
        except:
            desc = None

        cursor.execute('''SELECT a.msid, a.setkey, a.default_set, a.warning_low, 
                              a.caution_low, a.caution_high, a.warning_high, a.date, a.mlmenable, 
                              a.switchstate, a.mlimsw FROM limits AS a 
                              WHERE a.setkey = ? AND a.msid = ? AND a.datesec < ?
                              AND a.modversion = (SELECT MAX(b.modversion) FROM limits AS b
                              WHERE a.msid = b.msid and a.setkey = b.setkey and b.datesec < ?)''',
                       [setval, msid, DateTime(t1).secs, DateTime(t1).secs])
        b = cursor.fetchone()
        if not b:
            b = []

        cursor.execute('''SELECT a.msid, a.setkey, a.default_set, a.warning_low, 
                              a.caution_low, a.caution_high, a.warning_high, a.date, a.mlmenable, 
                              a.switchstate, a.mlimsw FROM limits AS a 
                              WHERE a.setkey = ? AND a.msid = ? AND a.datesec >= ? AND 
                              a.datesec <= ? AND a.modversion = (SELECT MAX(b.modversion) 
                              FROM limits AS b WHERE a.msid = b.msid AND a.setkey = b.setkey 
                              AND b.datesec >= ? AND b.datesec <= ?)''', [setval, msid, 
                              DateTime(t1).secs, DateTime(t2).secs, DateTime(t1).secs, 
                              DateTime(t2).secs])
        a = cursor.fetchone()

        changes[(msid, setval)] = {'before': b, 'after': a, 'description': desc}

    return changes


def write_report(thermal_msid_checks_file, t1, t2):
    """ Write Thermal Weekly Report

    :param thermal_msid_checks_file: Name of pickle file listing msids and related information
    :param t1: string containing start time in HOSC format
    :param t2: string containign stop time in HOSC format

    Note, in the past:
    thermal_msid_checks_file = 'thermalmsiddata.pkl'
    """
    thermdict, missing, notinarchive = pickle.load(open(thermal_msid_checks_file, 'r'))

    t1 = DateTime(t1).date
    t2 = DateTime(t2).date

    dayrange = (t1[:9] + '-' + t2[:9]).replace(':', '')

    power = get_average_tel_power(t1, t2)

    eclfile = 'ECLIPSE_HISTORY_2015.txt'
    ecltext = get_eclipse_text(eclfile, t1, t2)

    allviolations, missingmsids, checkedmsids = check_violations(thermdict, t1, t2)
    allviolations = add_violation_info(allviolations)

    # 3shtren and 4csdhav are not decommed correctly in the CXC archive
    if '3shtren' in allviolations.keys():
        _ = allviolations.pop('3shtren')

    if '4csdhav' in allviolations.keys():
        _ = allviolations.pop('4csdhav')

    limitchanges = check_limit_changes(t1, t2)

    html_limit_change_table = 'None'

    env = ja.Environment(loader=ja.FileSystemLoader('./templates'))

    template = env.get_template('thermal_weekly_template_.htm')
    webpage = template.render(startday=t1,
                              endday=t2,
                              dayrange=dayrange,
                              power=str('%5.1f' % power),
                              eclipse=ecltext,
                              violations=allviolations,
                              limitchanges=limitchanges)

    reportfilename = 'THERMAL_Weekly_' + dayrange + '.htm'
    outfile = file(reportfilename, 'w+')
    outfile.writelines(webpage)
    outfile.close()

    print('    Saved weekly report to {0}\n'.format(reportfilename))

    return reportfilename


def post_report(rootdir):

    dirs = os.listdir(rootdir)
    removed = [dirs.pop(dirs.index(d)) for d in dirs if len(d) != 15]
    for d in dirs:
        reportdir = pathjoin(rootdir, d)
        files = os.listdir(reportdir)
        t1str = d[:7]
        t2str = d[8:]
        t1 = DateTime(t1str + '.000000000').date
        t2 = DateTime(t2str + '.235959999').date

        thermalfile = 'THERMAL_Weekly_{}-{}.htm'.format(t1str, t2str)
        if thermalfile not in files:
            thermal_msid_checks_file = home + '/AXAFDATA/weekly_report_data/thermalmsiddata.pkl'
            newfile = write_report(thermal_msid_checks_file, t1, t2)
            copyfile(newfile, pathjoin(reportdir, newfile))
            copyfile('ThermalWeeklyReport.css', pathjoin(reportdir, 'ThermalWeeklyReport.css'))
            print 'copied files'


if __name__ == '__main__':

    #----------------------------------
    # Add and Parse Command Line Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--rootdir', default=None)
    args = vars(parser.parse_args())
    rootdir = args['rootdir']

    post_report(rootdir)



