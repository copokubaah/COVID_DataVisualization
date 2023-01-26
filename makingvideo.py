import numpy as np
import pandas as pd
import cv2 as cv

print('Libraries Imported')

cont = 'USA'
name = 'Deaths'
intpolnum = 90
where_to_start = 39

if cont == 'USA':
    compath = "C:/Users/copok/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"
    confirmed = pd.read_csv(compath)
    days = confirmed.columns[4:]
else:
    compath = "C:/Users/copok/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"
    confirmed = pd.read_csv(compath)
    days = confirmed.columns[4:]

fps = 60
out = cv.VideoWriter('D:/YoutubeChannel/CoronaVirus/{}_Bar{}_{}.avi'.format(name, cont, fps),
                     cv.VideoWriter_fourcc(*'DIVX'), fps, (1920, 1080))
months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN']
for dayidx, currday in enumerate(days):
    if dayidx >= where_to_start:
        if len(currday) == 6:
            month = months[int(currday[0]) - 1]
            day = currday[2]
        elif len(currday) == 7:
            month = months[int(currday[0]) - 1]
            day = currday[2:4]

        for num in np.arange(intpolnum):

            if (dayidx == where_to_start) and (num == 0):
                for _ in np.arange(fps*2):
                    figname = 'D:/YoutubeChannel/CoronaVirus/CoronaImages/Bar{}/{}{}-{}.png'.format(cont, month, day, num)
                    out.write(cv.imread(figname))
            elif (dayidx == (len(days)-1)) and (num == (intpolnum-1)):
                for _ in np.arange(fps*2):
                    figname = 'D:/YoutubeChannel/CoronaVirus/CoronaImages/Bar{}/{}{}-{}.png'.format(cont, month, day, num)
                    out.write(cv.imread(figname))
            else:
                figname = 'D:/YoutubeChannel/CoronaVirus/CoronaImages/Bar{}/{}{}-{}.png'.format(cont, month, day, num)
                out.write(cv.imread(figname))



        print('Day {} - {} Done'.format(dayidx, currday))

# close video
out.release()