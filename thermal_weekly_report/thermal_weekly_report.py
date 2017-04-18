
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


AXAFDATA = pathjoin(home, 'AXAFDATA')

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
    """Check a list of MSIDs for limit/expected state violations.
    
    :param thermdict: Dictionary of MSID information (MSID name, condition type, etc.)
    :param t1: String containing start date in HOSC format
    :param t2: String containgin stop date in HOSC format
    
    """
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
                
            if len(violations) > 0:
                allviolations[key] = process_violations(key, violations)
            
        except IndexError:
            print('{} not in DB'.format(key))
            missingmsids.append(key)

    return allviolations, missingmsids, checkedmsids

def handle_widerange_cases(key, t1, t2, greta_msid):
    """Handle special widerange MSIDs.
    
    :param key: Name of MSID as represented in Ska Engineering Archive
    :param t1: String containing start time in HOSC format
    :param t2: String containgin stop time in HOSC format
    :greta_msid: Name of MSID as represented in GRETA
    
    Note: Some MSID names differ between Ska and GRETA. Widerange MSIDs are one such case. For example
    OOBTHR35 is used for this measurement in both Ska and GRETA before this MSID was switched to 
    widerange read mode. Afterwards GRETA uses OOBTHR35_WIDE whereas Ska still uses OOBTHR35 for 
    continuity.
    """
    if DateTime(t2).secs <= DateTime('2014:342:16:30:00').secs:
        violations = pylimmon.check_limit_msid(key, t1, t2, greta_msid=key)
    elif DateTime(t1).secs >= DateTime('2014:342:16:33:00').secs:
        violations = pylimmon.check_limit_msid(key, t1, t2, greta_msid=greta_msid)
    else:
        t2_a = np.min((DateTime(t2).secs, DateTime('2014:342:16:30:00').secs))
        violations = pylimmon.check_limit_msid(key, t1, t2_a, greta_msid=key)
        t1_b = np.min((DateTime(t2).secs, DateTime('2014:342:16:33:00').secs))
        violations_b = pylimmon.check_limit_msid(key, t1_b, t2, greta_msid=greta_msid)

        violations.extend(violations_b)
        
    return violations

def process_violations(msid, violations):
    """Add contextual information for any limit/expected state violations.
    
    :param msid: Current mnemonic
    :param violations: List of individual violations (list of tuples)
    
    """
    data = fetch.Msid(msid, violations[0][0][0], violations[0][0][-1], stat='5min')
    try:
        desc = data.tdb.technical_name
    except:
        desc = 'No Description in TDB'
        
    violation_dict = {}
    for v in violations:
        limtype = v[-1]
        if 'high' in limtype.lower():
            if limtype not in violation_dict.keys():
                violation_dict.update({limtype:{'starttime':v[0][0], 'stoptime':v[0][-1], 'num_excursions':1,
                                                'extrema':np.max(v[1]), 'limit':v[2][0], 'setid':v[3][0],
                                                'duration':v[0][-1] - v[0][0]}})
            else:
                violation_dict[limtype]['extrema'] = np.max((np.max(v[1]), violation_dict[limtype]['extrema']))
                violation_dict[limtype]['starttime'] = np.min((v[0][0], violation_dict[limtype]['starttime']))
                violation_dict[limtype]['stoptime'] = np.max((v[0][0], violation_dict[limtype]['stoptime']))
                violation_dict[limtype]['num_excursions'] = violation_dict[limtype]['num_excursions'] + 1
                violation_dict[limtype]['duration'] = violation_dict[limtype]['duration'] + v[0][-1] - v[0][0]

        elif 'low' in limtype.lower():
            if limtype not in violation_dict.keys():
                violation_dict.update({limtype:{'starttime':v[0][0], 'stoptime':v[0][-1], 'num_excursions':1,
                                                'extrema':np.min(v[1]), 'limit':v[2][0], 'setid':v[3][0],
                                                'duration':v[0][-1] - v[0][0]}})
            else:
                violation_dict[limtype]['extrema'] = np.min((np.min(v[1]), violation_dict[limtype]['extrema']))
                violation_dict[limtype]['starttime'] = np.min((v[0][0], violation_dict[limtype]['starttime']))
                violation_dict[limtype]['stoptime'] = np.max((v[0][0], violation_dict[limtype]['stoptime']))
                violation_dict[limtype]['num_excursions'] = violation_dict[limtype]['num_excursions'] + 1
                violation_dict[limtype]['duration'] = violation_dict[limtype]['duration'] + v[0][-1] - v[0][0]

        elif 'state' in limtype.lower():
            if limtype not in violation_dict.keys():
                violation_dict.update({limtype:{'starttime':v[0][0], 'stoptime':v[0][-1], 'num_excursions':1,
                                                'extrema':v[1][0], 'limit':v[2][0], 'setid':v[3][0],
                                                'duration':v[0][-1] - v[0][0]}})
            else:
                violation_dict[limtype]['starttime'] = np.min((v[0][0], violation_dict[limtype]['starttime']))
                violation_dict[limtype]['stoptime'] = np.max((v[0][0], violation_dict[limtype]['stoptime']))
                violation_dict[limtype]['num_excursions'] = violation_dict[limtype]['num_excursions'] + 1
                violation_dict[limtype]['duration'] = violation_dict[limtype]['duration'] + v[0][-1] - v[0][0]

                
    for limittype in ['warning_low', 'caution_low', 'caution_high', 'warning_high', 'state']:
        if limittype in violation_dict.keys():
            violation_dict[limittype]['duration'] = violation_dict[limittype]['duration'] / 3600.
            violation_dict[limittype]['description'] = desc
            violation_dict[limittype]['startdate'] = DateTime(violation_dict[limittype]['starttime']).date
            violation_dict[limittype]['stopdate'] = DateTime(violation_dict[limittype]['stoptime']).date
            
    return violation_dict



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

    eclfile = pathjoin(AXAFDATA, 'ECLIPSE.txt')
    ecltext = get_eclipse_text(eclfile, t1, t2)

    allviolations, missingmsids, checkedmsids = check_violations(thermdict, t1, t2)

    # 3shtren and 4csdhav are not decommed correctly in the CXC archive
    if '3shtren' in allviolations.keys():
        _ = allviolations.pop('3shtren')

    if '4csdhav' in allviolations.keys():
        _ = allviolations.pop('4csdhav')

    limitchanges = check_limit_changes(t1, t2)

    html_limit_change_table = 'None'

    env = ja.Environment(loader=ja.FileSystemLoader(home + '/AXAFLIB/thermal_weekly_report/thermal_weekly_report/templates'))

    template = env.get_template('thermal_weekly_template.htm')
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
            copyfile(home + '/AXAFLIB/thermal_weekly_report/thermal_weekly_report/ThermalWeeklyReport.css', pathjoin(reportdir, 'ThermalWeeklyReport.css'))
            print 'copied files'


if __name__ == '__main__':

    #----------------------------------
    # Add and Parse Command Line Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--rootdir', default=None)
    args = vars(parser.parse_args())
    rootdir = args['rootdir']

    post_report(rootdir)



