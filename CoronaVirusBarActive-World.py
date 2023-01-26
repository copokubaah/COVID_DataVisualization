## CREATING ANIMATION FOR MAP AFRICA

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as color
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import cv2 as cv
import geopandas as gpd
import matplotlib.pyplot as plt
from scipy import ndimage
import matplotlib.patheffects as path_effects
import matplotlib as mpl
import mapclassify
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import gc
import concurrent.futures
from unicodedata import normalize

print('Libraries Imported')


# ## Create Functions

def cleandata(data):
    p1 = data.groupby('Country/Region').first()
    p1 = p1[['Lat', 'Long']]
    p2 = data.groupby('Country/Region').sum()
    return p1.join(p2.iloc[:, 2:]).reset_index()


def get_image(name):
    path = "D:/YoutubeChannel/CoronaVirus/{}".format(name)
    im = plt.imread(path)
    im = (im * 255).astype(np.uint8)
    return im

def get_image2(name):
    path = "D:/YoutubeChannel/CoronaVirus/{}".format(name)
    im = plt.imread(path)
    # im = (im * 255).astype(np.uint8)
    return im


def int2str(text):
    if len(text) <= 3:
        text1 = text
    elif (len(text) > 3) and (len(text) < 7):
        text1 = ','.join([text[:-3], text[-3:]])
    elif (len(text) > 6) and (len(text) < 10):
        text1 = ','.join([text[:-6], text[-6:-3], text[-3:]])
    elif len(text) > 9:
        text1 = ','.join([text[:-9], text[-9:-6], text[-6:-3], text[-3:]])
    return text1

def offset_image(xcoord, ycoord, name, ax, fig, zoomscalar, xres=600, yres=400):
    img1 = get_image(name)
    img = cv.resize(img1, (xres, yres), interpolation=cv.INTER_AREA)
    im = OffsetImage(img, zoom=zoomscalar, transform=ax.transAxes)
    im.image.axes = ax
    # rend = fig.canvas.get_renderer()
    # xadd = im.get_extent(rend)[0]
    ab = AnnotationBbox(im, (xcoord, ycoord), xybox=(0., 0.), frameon=False,
                        xycoords='data', boxcoords="offset points", pad=0)
    ax.add_artist(ab)


def offset_image2(xcoord, ycoord, img1, ax, zoomscalar, xres=600, yres=400):
    img = cv.resize(img1, (xres, yres), interpolation=cv.INTER_AREA)
    im = OffsetImage(img, zoom=zoomscalar, transform=ax.transAxes)
    im.image.axes = ax
    ab = AnnotationBbox(im, (xcoord, ycoord), xybox=(0., 0.), frameon=False,
                        xycoords='data', boxcoords="offset points", pad=0)
    ax.add_artist(ab)


def interpol_values(x0, x1, intpolnum):
    rnum = len(x0)
    xall = np.zeros((rnum, intpolnum))
    for num in np.arange(rnum):
        xall[num, :] = np.linspace(x0[num], x1[num], intpolnum, endpoint=False)
    return xall.astype(int)


def all_intpol(days, where_to_start, data, intpolnum):
    daysum = len(days)
    intdata = np.empty((data.shape[0], 1))
    for dayidx, currday in enumerate(days):
        if dayidx >= where_to_start:
            if dayidx == daysum - 1:
                nextday = currday
            else:
                nextday = days[dayidx + 1]
            intdata = np.hstack((intdata, interpol_values(data[currday].values, data[nextday].values, intpolnum)))
    return intdata[:, 1:]


def daylist(days, where_to_start, intpolnum, months=['JAN', 'FEB', 'MAR', 'APR', 'MAY']):
    day = []
    month = []
    for dayidx, currday in enumerate(days):
        if dayidx >= where_to_start:
            if len(currday) == 6:
                for i in np.arange(intpolnum):
                    month.append(months[int(currday[0]) - 1])
                    day.append(currday[2])
            elif len(currday) == 7:
                for i in np.arange(intpolnum):
                    month.append(months[int(currday[0]) - 1])
                    day.append(currday[2:4])
    return day, month


def smooth_transition(data, numt, topnum):
    ind = np.argsort(data, axis=0)
    pos = np.zeros((ind.shape[0], ind.shape[1]))
    for i in np.arange(ind.shape[0]):
        for j in np.arange(ind.shape[1]):
            pos[i, j] = np.nonzero(ind[:, j] == i)[0][0]

    smoth = np.copy(pos)
    cutval = pos.shape[0]-21
    for i in np.arange(pos.shape[0]):
        d = np.diff(pos[i, :])
        dn = np.nonzero(d != 0)[0]

        for j in dn:
            if pos[i,j] < (pos.shape[0]-topnum) and pos[i,j+1] >= (pos.shape[0]-topnum):
                smoth[i, j:(j + numt)] = np.linspace(cutval, pos[i, j + 1], numt)
            elif pos[i,j+1] < (pos.shape[0]-topnum) and pos[i,j] >= (pos.shape[0]-topnum):
                smoth[i, j:(j + numt)] = np.linspace(pos[i, j], cutval, numt)
            elif pos[i,j+1] < (pos.shape[0]-topnum) and pos[i,j] < (pos.shape[0]-topnum):
                smoth[i, j:(j + numt)] = np.linspace(cutval, cutval, numt)
            elif pos[i, j + 1] >= (pos.shape[0] - topnum) and pos[i, j] >= (pos.shape[0] - topnum):
                smoth[i, j:(j + numt)] = np.linspace(pos[i, j], pos[i, j + 1], numt)
    return ind, pos, smoth

def changestate(confirmed, death, rec, chdict):
    for t in np.arange(len(chdict)):
        confirmed.loc[confirmed['Country/Region'] == chdict[t][0], 'Country/Region'] = chdict[t][1]
        death.loc[death['Country/Region'] == chdict[t][0], 'Country/Region'] = chdict[t][1]
        rec.loc[rec['Country/Region'] == chdict[t][0], 'Country/Region'] = chdict[t][1]
    return confirmed, death, rec

def changestate2(confirmed, chdict):
    for t in np.arange(len(chdict)):
        confirmed.loc[confirmed['Country/Region'] == chdict[t][0], 'Country/Region'] = chdict[t][1]
    return confirmed


# ## Import and cleaning data

# compath = "D:/YoutubeChannel/CoronaVirus/time_series_covid19_confirmed_global.csv"
# dpath = "D:/YoutubeChannel/CoronaVirus/time_series_covid19_deaths_global.csv"
# rpath = "D:/YoutubeChannel/CoronaVirus/time_series_covid19_recovered_global.csv"
compath = "C:/Users/copok/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
dpath = "C:/Users/copok/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
rpath = "C:/Users/copok/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"
confirmed = pd.read_csv(compath)
death = pd.read_csv(dpath)
rec = pd.read_csv(rpath)
print('The shapes of confirmed, death and recovered are {}, {} and {}'.format(confirmed.shape, death.shape, rec.shape))

confirmed.loc[confirmed['Province/State'] == 'Greenland', 'Country/Region'] = 'Greenland'
death.loc[death['Province/State'] == 'Greenland', 'Country/Region'] = 'Greenland'
rec.loc[rec['Province/State'] == 'Greenland', 'Country/Region'] = 'Greenland'
chdict = [['Korea, South','South Korea'], ['Taiwan*','Taiwan'], ['Burma', 'Myanmar']]
confirmed, death, rec = changestate(confirmed, death, rec, chdict)

conf = cleandata(confirmed)
deat = cleandata(death)
recd = cleandata(rec)
conf = conf.fillna(0)
deat = deat.fillna(0)
rec = rec.fillna(0)
print('The shapes of confirmed, death and recovered are {}, {} and {}'.format(conf.shape, deat.shape, recd.shape))

popath = "D:/YouTubeChannel/Technology/population_total.csv"
population = pd.read_csv(popath)
print('The shapes of confirmed is {}'.format(confirmed.shape))

chdict = [['Congo, Dem. Rep.', 'Congo (Kinshasa)'], ['Congo, Rep.', 'Congo (Brazzaville'], ['Czech Republic', 'Czechia'],
          ['Kyrgyz Republic', 'Kyrgyzstan'], ['Lao', 'Laos'], ['Micronesia, Fed. Sts.', 'Micronesia'],
          ['Slovak Republic', 'Slovakia'], ['St. Lucia', 'Saint Lucia'], ['Swaziland', 'Eswatini'],
          ['St. Vincent and the Grenadines', 'Saint Vincent and the Grenadines'], ['United States', 'US']]
population = population.rename(columns={'Country': 'Country/Region'})
population = changestate2(population, chdict)

uscode = pd.read_csv("D:/YoutubeChannel/CoronaVirus/worldcodes2.csv")
uscode = uscode.rename(columns={'Country': 'Country/Region'})
popuscode = pd.merge(population[['Country/Region', '2020']], uscode, on='Country/Region', how='inner')


# download us codes and choose the continent


popuscode.loc[popuscode['Country/Region'] == 'Namibia', 'code'] = 'NA'
sort_conf = pd.merge(popuscode[['code', 'Continent_Name', 'Country/Region', '2020']], conf, on='Country/Region', how='inner')
sort_death = pd.merge(popuscode[['code', 'Continent_Name', 'Country/Region', '2020']], deat, on='Country/Region', how='inner')
sort_rec = pd.merge(popuscode[['code', 'Continent_Name', 'Country/Region', '2020']], recd, on='Country/Region', how='inner')

# sort_death = sort_death.loc[sort_conf[sort_conf.columns[-1]] > 1000, :].reset_index(drop=True)
# sort_rec = sort_rec.loc[sort_conf[sort_conf.columns[-1]] > 1000, :].reset_index(drop=True)
# sort_conf = sort_conf.loc[sort_conf[sort_conf.columns[-1]] > 1000, :].reset_index(drop=True)

print('The shapes of confirmed, death and recovered are {}, {} and {}'.format(sort_conf.shape,
                                                                              sort_death.shape, sort_rec.shape))
conlist = list(sort_conf.iloc[:, 6:].max())
max_cases = max(conlist) + 100
print('The maximum cases is {}'.format(max_cases))

# ## Import Accessory Data such as Flags, Images and world map
# download situation reports and flags
report = pd.read_excel("D:/YoutubeChannel/CoronaVirus/SituationReports.xlsx")
flagpath = "D:/YoutubeChannel/CoronaVirus/worldflags.csv"
flagname = pd.read_csv(flagpath)

# get coronavirus image
corimg = get_image('coronavirus.png')
bgimg = get_image2('worldmap/w39.jpg')
# cutlum, replum = 20,10
# bgimg[bgimg<cutlum] = bgimg[bgimg<cutlum] + replum
bgimg = bgimg + 10
bgimg = cv.resize(bgimg, (1920, 1080), interpolation=cv.INTER_AREA)




# ## Create Parent Class
# #### The parent class will contain class variables comprising of all cleaned data, data accessorys and plot variables

class Create_Animation():
    plot_dict = None
    acc_dict = None

    def __init__(self, data_dict, datatype):
        self.data_dict = data_dict
        self.datatype = datatype

    @classmethod
    def set_variables(cls, plot_dict):
        cls.plot_dict = plot_dict

    def intpoldata(self, intpolnum, topnum, numt):
        self.intpolnum = intpolnum
        self.numt = numt
        self.topnum = topnum
        self.per = 1000000
        self.days = self.data_dict['sort_conf'].columns[6:]
        self.conf_all = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['sort_conf'],
                                   self.intpolnum)
        self.death_all = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['sort_death'],
                                    self.intpolnum)
        self.rec_all = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['sort_rec'],
                                  self.intpolnum)
        self.act_all = self.conf_all - (self.death_all + self.rec_all)
        self.conf_w = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['conf'], self.intpolnum)
        self.death_w = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['deat'], self.intpolnum)
        self.rec_w = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['recd'], self.intpolnum)
        self.day, self.month = daylist(self.days, self.plot_dict['where_to_start'], self.intpolnum)
        self.act_w = self.conf_w - (self.death_w + self.rec_w)


        if self.datatype == 'Total':
            self.data1 = self.conf_all
            self.data2 = self.data_dict['sort_conf'][self.data_dict['sort_conf'].columns[:3]].join(
                pd.DataFrame(self.conf_all, columns=self.day))
            self.data_w = self.conf_w
        elif self.datatype == 'Active':
            self.data1 = self.act_all
            self.data2 = self.data_dict['sort_conf'][self.data_dict['sort_conf'].columns[:4]].join(
                pd.DataFrame(self.act_all, columns=self.day))
            self.data_w = self.act_w
        elif self.datatype == 'Deaths':
            self.data1 = self.death_all
            self.data2 = self.data_dict['sort_conf'][self.data_dict['sort_conf'].columns[:4]].join(
                pd.DataFrame(self.death_all, columns=self.day))
            self.data_w = self.death_w
        elif self.datatype == 'Recovered':
            self.data1 = self.rec_all
            self.data2 = self.data_dict['sort_conf'][self.data_dict['sort_conf'].columns[:4]].join(
                pd.DataFrame(self.rec_all, columns=self.day))
            self.data_w = self.rec_w

        self.pop = self.data2.iloc[:, 3].values

        self.africa = self.data2[self.data2.Continent_Name == 'Africa']
        afcon = list(self.africa['Country/Region'])
        afpop = self.africa.iloc[:, 3].values
        self.africa = self.africa.iloc[:, 4:].values

        self.asia = self.data2[self.data2.Continent_Name == 'Asia']
        ascon = list(self.asia['Country/Region'])
        aspop = self.asia.iloc[:, 3].values
        self.asia = self.asia.iloc[:, 4:].values
        self.america = self.data2.loc[((self.data2['Continent_Name'] == 'North America') | (
                    self.data2['Continent_Name'] == 'South America')), :]
        amcon = list(self.america['Country/Region'])
        ampop = self.america.iloc[:, 3].values
        self.america = self.america.iloc[:, 4:].values
        self.europe = self.data2[self.data2.Continent_Name == 'Europe']
        eucon = list(self.europe['Country/Region'])
        eupop = self.europe.iloc[:, 3].values
        self.europe = self.europe.iloc[:, 4:].values
        self.oceania = self.data2[self.data2.Continent_Name == 'Oceania']
        occon = list(self.oceania['Country/Region'])
        ocpop = self.oceania.iloc[:, 3].values
        self.oceania = self.oceania.iloc[:, 4:].values

        # conafrica= self.conf2[self.conf2.Continent_Name == 'Africa'].iloc[:, 3:].values
        # conasia = self.conf2[self.conf2.Continent_Name == 'Asia'].iloc[:, 3:].values
        # conamerica = self.conf2.loc[((self.conf2['Continent_Name'] == 'North America') | (
        #             self.conf2['Continent_Name'] == 'South America')), :].iloc[:, 3:].values
        # coneurope = self.conf2[self.conf2.Continent_Name == 'Europe'].iloc[:, 3:].values
        # conoceania = self.conf2[self.conf2.Continent_Name == 'Oceania'].iloc[:, 3:].values

        self.wsum = np.sum(self.data1, axis=0)
        afsum = np.sum(self.africa, axis=0)
        assum = np.sum(self.asia, axis=0)
        amsum = np.sum(self.america, axis=0)
        eusum = np.sum(self.europe, axis=0)
        ocsum = np.sum(self.oceania, axis=0)

        afall = np.true_divide(afsum, np.sum(afpop)) * self.per
        asall = np.true_divide(assum, np.sum(aspop)) * self.per
        amall = np.true_divide(amsum, np.sum(ampop)) * self.per
        euall = np.true_divide(eusum, np.sum(eupop)) * self.per
        ocall = np.true_divide(ocsum, np.sum(ocpop)) * self.per
        self.wall = np.true_divide(self.wsum, np.sum(self.pop)) * self.per

        self.contorder = np.vstack((afall.reshape(1, -1), amall.reshape(1, -1),
                                    asall.reshape(1, -1), euall.reshape(1, -1)))
        self.contsum = np.vstack((afsum.reshape(1, -1), amsum.reshape(1, -1),
                                    assum.reshape(1, -1), eusum.reshape(1, -1)))

        if self.plot_dict['which_cont'] == 'Africa':
            self.data = np.true_divide(self.africa, afpop.reshape(-1, 1)) * self.per
            self.cntryname = afcon
        elif self.plot_dict['which_cont'] == 'Asia':
            self.data = np.true_divide(self.asia, aspop.reshape(-1, 1)) * self.per
            self.cntryname = ascon
        elif self.plot_dict['which_cont'] == 'America':
            self.data = np.true_divide(self.america, ampop.reshape(-1, 1)) * self.per
            self.cntryname = amcon
        elif self.plot_dict['which_cont'] == 'Europe':
            self.data = np.true_divide(self.europe, eupop.reshape(-1, 1)) * self.per
            self.cntryname = eucon
        elif self.plot_dict['which_cont'] == 'World':
            self.data = np.nan_to_num(np.true_divide(self.data1, self.pop.reshape(-1, 1)) * self.per)
            self.cntryname = list(k.data_dict['sort_conf']['Country/Region'])

        self.ind, self.pos, self.smoth = smooth_transition(self.data, self.numt, self.topnum)
        self.indcont, self.poscont, self.smothcont = smooth_transition(self.contorder, self.numt, 4)

    def getbardata(self, num):

        self.num = num
        self.numname = self.num % self.intpolnum
        self.topidx = np.flip(self.ind[:, self.num])[:self.topnum]
        self.ynames = [self.cntryname[i] for i in self.topidx]
        self.y = self.smoth.shape[0] - self.smoth[self.topidx, self.num]
        self.xdata = self.data[self.topidx, self.num]
        # self.xact = self.act_all[self.topidx, self.num]
        # self.xconf = self.conf_all[self.topidx, self.num]
        # self.xdeath = self.death_all[self.topidx, self.num]
        # self.xrec = self.rec_all[self.topidx, self.num]

        # determine if topidx contains zeros
        testconf = self.data[self.topidx, self.num]
        self.numsigdx = len(np.nonzero(testconf)[0])

    def plotbargraph(self, title, rect, bgimg):

        # set figure and axis
        fig = plt.figure(figsize=(self.plot_dict['myres_w'] / self.plot_dict['mydpi'],
                                  self.plot_dict['myres_h'] / self.plot_dict['mydpi']),
                         dpi=self.plot_dict['mydpi'], facecolor=self.plot_dict['bkgcolor'])

        self.fig = fig
        self.axb = self.fig.add_axes([0, 0, 1, 1])
        self.axb.imshow(bgimg)
        self.axb.set_axis_off()

        self.ax = self.fig.add_axes(rect)
        self.h = self.ax.barh(self.y, self.xdata, height=0.7)
        self.maxdata = np.max(self.xdata)
        xlimit = self.maxdata * 1.05
        scal = np.array([i * (10**(len(str(int(xlimit)))-2)) for i in [1, 2, 5, 10, 20, 50]])
        ngrid = np.array([xlimit//i for i in scal])
        stepgrid = np.min(scal[ngrid < 10])

        contcolor = []
        for i, cntry in enumerate(self.ynames):
            if i < self.numsigdx:
                # print(cntry)
                cncode = list(self.data_dict['sort_conf'].loc[self.data_dict['sort_conf']['Country/Region'] == cntry, 'code'])[0].lower()
                contcolor.append(list(flagname.loc[flagname['code'] == cncode, 'color'])[0])
                self.h[i].set_color(contcolor[i])


        # set axis properties
        self.ax.set_ylim(0, self.topnum + 1)
        self.ax.set_yticks(self.y[:self.numsigdx])
        self.ax.set_yticklabels(self.ynames[:self.numsigdx],
                                fontdict={'fontname': self.plot_dict['myfont'], 'fontweight': 'semibold', 'fontsize': 5,
                                          'va': 'center'})
        for ytick, color in zip(self.ax.get_yticklabels(), ['w']*20): #contcolor
            ytick.set_color(color)
        self.ax.tick_params(axis='both', which='both', length=0, pad=2)
        xtickrange = np.round((np.arange(0, self.maxdata, stepgrid))).astype('int')
        self.ax.set_xticks(xtickrange)
        self.ax.set_xlim(0, 1.2 * xlimit)
        self.ax.set_xticklabels(xtickrange,
                                fontdict={'fontname': self.plot_dict['myfont'], 'fontsize': 4, 'fontweight': 3},
                                color=(0.5, 0.5, 0.5, 0.8))
        plt.grid(axis='x', c='gray', ls='--', lw=0.6, alpha=0.3)
        self.ax.set_title(title, x=-0.05, fontdict={'fontsize': 15, 'fontname': self.plot_dict['myfont'], 'fontweight': 'bold',
                                           'ha': 'left'},
                          color='deepskyblue', path_effects=[
                path_effects.withSimplePatchShadow(offset=(2, -2), shadow_rgbFace='k', alpha=0.3)])


    def add_items(self, zoomscalar = 0.018):

        # insert flag and annotate the number of cases at the right of it
        self.zoomscalar = zoomscalar
        imagedist = self.maxdata / 38
        for i, cntry in enumerate(self.ynames):
            if i < self.numsigdx:
                cncode = list(self.data_dict['sort_conf'].loc[self.data_dict['sort_conf']['Country/Region'] == cntry, 'code'])[0].lower()
                c = f'worldflag/{cncode}.png'
                offset_image(self.xdata[i] + imagedist, self.y[i], c, self.ax, self.fig, self.zoomscalar)

        for py, px in enumerate(self.xdata):
            text1 = np.round(px, 2)
            contname = list(self.data_dict['sort_conf'].loc[
                                self.data_dict['sort_conf']['Country/Region'] == self.ynames[py], 'Continent_Name'])[
                0].upper()
            if contname == 'NORTH AMERICA':
                contname = 'N. AMERICA'
            elif contname == 'SOUTH AMERICA':
                contname = 'S. AMERICA'
            if py < self.numsigdx:
                self.ax.annotate(text1, xy=(px + imagedist + self.maxdata / 38, self.y[py]),
                                 fontname=self.plot_dict['myfont'],
                                 fontsize=4.2, fontweight='semibold', fontstyle='normal', va='center', ha='left',
                                 color='w')
                if self.plot_dict['which_cont'] == 'World':
                    if (len(contname) > 9):
                        if ((px / self.maxdata) > 0.09):
                            self.ax.annotate(contname, xy=(px - self.maxdata / 250, self.y[py]), color='w',
                                             fontweight='semibold',
                                             fontname=self.plot_dict['myfont'], fontsize=4.2, fontstyle='normal',
                                             va='center', ha='right')
                    else:
                        if ((px / self.maxdata) > 0.06):
                            self.ax.annotate(contname, xy=(px - self.maxdata / 250, self.y[py]), color='w',
                                             fontweight='semibold',
                                             fontname=self.plot_dict['myfont'], fontsize=4.2, fontstyle='normal',
                                             va='center', ha='right')


    def add_text(self, sidetop, rect1=[0.84, 0, 0.15, 0.85]):
        self.ax1 = self.fig.add_axes(rect1)

        headcolor = 'w'
        xs = 0.4
        hour = int(int(self.num % self.intpolnum) // (self.intpolnum / 24))
        minute = int(int(self.num % self.intpolnum) // (self.intpolnum / 60))  #int((self.num % (self.intpolnum / 24)) * 60 / (self.intpolnum / 24))
        self.ax1.text(xs, 0.83, '{}:{}'.format(hour, minute), transform=self.ax1.transAxes,fontweight='semibold',
                     verticalalignment='center', horizontalalignment='center', color='r', fontsize=10,
                     fontname=self.plot_dict['myfont'])
        self.ax1.text(xs, 0.88, '{} {}'.format(self.day[self.num], self.month[self.num]), transform = self.ax1.transAxes,fontweight='semibold',
                     verticalalignment='center', horizontalalignment='center', color='w', fontsize=18, fontname=self.plot_dict['myfont'])
        self.ax1.text(xs, 0.74, f'Active Cases\nby Continent', transform=self.ax1.transAxes, fontweight='semibold',
                      verticalalignment='center', horizontalalignment='center', color=headcolor, fontsize=7,
                      fontname=self.plot_dict['myfont'])
        self.ax1.text(xs, 0.685, f'(% of World Total)', transform=self.ax1.transAxes,
                      verticalalignment='center', horizontalalignment='center', color=headcolor, fontsize=5,
                      fontweight='semibold',
                      fontname=self.plot_dict['myfont'])

        ystp = 0.05
        ysd = 0.12
        bigfont = 8
        smallfont = 6
        contname = ['Africa', 'America', 'Asia', 'Europe']
        contcolor = ['green', 'royalblue', 'firebrick', 'gold']
        for n, idx, i in zip(np.arange(len(contname)), self.indcont[:, self.num], self.smothcont[:, self.num]):
            ys = 0.64 - (len(contname)-1-n) * ysd

            tnum = np.round(self.contorder[idx, self.num], 2)
            tsum = np.round(self.contsum[idx, self.num])
            text = str(int(tsum))
            text1 = int2str(text)
            cntry = contname[idx]

            self.ax1.text(xs, ys, '{}'.format(cntry.upper()), transform=self.ax1.transAxes, fontweight='semibold',
                          va='center', ha='center', color=contcolor[idx], fontsize=bigfont,
                          fontname=self.plot_dict['myfont'])
            self.ax1.text(xs, ys - ystp, '{}\n{:2.2f} per 1M pop'.format(text1, tnum), fontweight='semibold',
                          transform=self.ax1.transAxes, va='center', ha='center', color=contcolor[idx],
                          fontsize=smallfont, fontname=self.plot_dict['myfont'])

        ys = 0.64 - len(contname) * ysd
        tworld = self.wall[self.num]
        tworldsum = self.wsum[self.num]
        text = str(int(tworldsum))
        text1 = int2str(text)

        self.ax1.text(xs, ys, 'WORLD', transform=self.ax1.transAxes,
                      va='center', ha='center', color='w', fontweight='semibold', fontsize=bigfont + 2,
                      fontname=self.plot_dict['myfont'])
        self.ax1.text(xs, ys - ystp*1.2, '{}\n{:2.2f} per 1M pop'.format(text1, tworld), transform=self.ax1.transAxes, va='center', ha='center',
                      color='w', fontsize=smallfont + 2, fontweight='semibold', fontname=self.plot_dict['myfont'])

        self.ax1.set_axis_off()

    def add_report(self, report):

        # insert coronavirus image
        self.ax2 = self.fig.add_axes([0.25, 0.15, 0.2, 0.2])
        degrees = np.linspace(0, 360, 5*self.intpolnum)
        rotated_img = ndimage.rotate(corimg, degrees[int(self.num % (5*intpolnum))])
        offset_image2(0.5, 0.5, corimg, self.ax2, 0.12, 500, 500)
        self.ax2.set_axis_off()

        xst = 0.37
        yst = 0.25
        edcolor = 'black'
        talpha = 0
        tcolor = (0.9, 0.9, 0.9, 0.7)
        textlim = 40
        dayidx = int(self.num // self.intpolnum) + self.plot_dict['where_to_start']
        if self.plot_dict['which_cont'] == 'Africa':
            which_report = 2
        else:
            which_report = 1

        rpnum = int(type(report.iloc[dayidx, which_report]) == str)
        if rpnum == 1:
            if self.numname > int(self.intpolnum * 0.1) and self.numname < int(self.intpolnum * 0.9):
                text1 = normalize('NFKD', report.iloc[dayidx, which_report])
                tx = text1.split()
                cumtx = np.cumsum(np.array([len(i)+1 for i in tx]))
                self.newtxt = ''
                n = 1

                while len(text1) > len(self.newtxt):
                    lastindx = np.int(np.nonzero(cumtx < textlim*n)[0][-1])
                    if n == 1:
                        self.newtxt = " ".join(tx[: lastindx+1])
                    else:
                        self.newtxt = self.newtxt + '\n' + " ".join(tx[previndx+1: lastindx+1])

                    previndx = lastindx
                    n += 1

                self.ax.text(xst, yst, '{}'.format(self.newtxt), transform=self.ax.transAxes,
                             verticalalignment='center', horizontalalignment='left', color=tcolor, fontsize=7,
                             fontname=self.plot_dict['myfont'],fontweight='semibold',
                             bbox=dict(fc='w', ec=edcolor, pad=5, alpha=talpha))


    def set_finalaxis(self):
        self.ax.invert_yaxis()
        self.ax.set_frame_on(False)


    def saveplot(self):
        # write to video
        figname = 'D:/YoutubeChannel/CoronaVirus/CoronaImages/Bar{}/{}{}-{}.png'.format(self.plot_dict['which_cont'], self.month[self.num],
                                                                                         self.day[self.num], self.numname)
        self.fig.savefig(figname, dpi=self.fig.get_dpi(), facecolor=self.fig.get_facecolor())
        plt.cla()
        self.fig.clf()
        plt.close('all')
        gc.collect()


def plotgraph(num):

    k.getbardata(num)
    if which_cont == 'World':
         title = f'Top 20 Countries by {datatype} COVID-19 Cases'
    else:
        title = f'Top 20 Countries in {which_cont} by {datatype} COVID-19 Cases'
    rect = [0.15, 0.05, 0.8, 0.85]
    k.plotbargraph(title, rect, bgimg)
    zoomscalar = 0.02
    k.add_items(zoomscalar)
    sidetop = 3
    k.add_text(sidetop)
    k.add_report(report)
    k.set_finalaxis()
    k.saveplot()

    return f'done with {num}'

################################################################################################################
where_to_start = 39
startnow = 0
data_dict = dict(sort_conf=sort_conf, sort_death=sort_death, sort_rec=sort_rec, conf=conf, deat=deat, recd=recd)
myfont = 'Consolas'
bkgcolor = 'ghostwhite'
which_cont = 'Africa'
plot_dict = dict(where_to_start=where_to_start, myfont=myfont, bkgcolor=bkgcolor, myres_w=1920,
                 myres_h=1080, mydpi=300, which_cont=which_cont)

datatype='Deaths'
k = Create_Animation(data_dict, datatype)
Create_Animation.set_variables(plot_dict)
intpolnum = 210
topnum = 20
numt = 60
k.intpoldata(intpolnum, topnum, numt)
totalnum = intpolnum * (len(sort_conf.columns[6:]) - where_to_start)
startnum = 13440   #intpolnum * startnow

num = totalnum-100
plotgraph(num)

#################
#
# if __name__ == '__main__':
#
#      with concurrent.futures.ProcessPoolExecutor() as executor:
#         num = range(startnum, totalnum)
#         results = executor.map(plotgraph, num)
#
#         for result in results:
#             print(result)




