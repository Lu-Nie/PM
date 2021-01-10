from datetime import datetime
 
 
def d_to_jd(time):
    fmt = '%Y.%m.%d'
    dt = datetime.strptime(time, fmt)
    tt = dt.timetuple()
    return tt.tm_year * 1000 + tt.tm_yday
 
 
def jd_to_time(time):
    dt = datetime.strptime(time, '%Y%j').date()
    fmt = '%Y.%m.%d'
    return dt.strftime(fmt)
 
 
if '__main__' == __name__:
    flag = input('Please input the date transform type: 1 for Julian day to date; 2 for date to Julian day:')
    if 1 == int(flag):
        time = input('Please input the date (YY.MM.DD):')
        print(d_to_jd(time))
    elif 2 == int(flag):
        time = input('Please input the Julian day:')
        print(jd_to_time(time))