## CREATING ANIMATION FOR MAP AFRICA

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import cv2 as cv
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.patheffects as path_effects
from scipy import ndimage
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
    return xall


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


def daylist(days, where_to_start, intpolnum, months=['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN']):
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
            self.data2 = self.data_dict['sort_conf'][self.data_dict['sort_conf'].columns[:4]].join(
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

        if self.plot_dict['which_cont'] == 'Africa':
            self.data = np.true_divide(self.africa, afpop.reshape(-1, 1)) * self.per
            self.rdata = self.africa
            self.cntryname = afcon
        elif self.plot_dict['which_cont'] == 'Asia':
            self.data = np.true_divide(self.asia, aspop.reshape(-1, 1)) * self.per
            self.rdata = self.asia
            self.cntryname = ascon
        elif self.plot_dict['which_cont'] == 'America':
            self.data = np.true_divide(self.america, ampop.reshape(-1, 1)) * self.per
            self.rdata = self.america
            self.cntryname = amcon
        elif self.plot_dict['which_cont'] == 'Europe':
            self.data = np.true_divide(self.europe, eupop.reshape(-1, 1)) * self.per
            self.rdata = self.europe
            self.cntryname = eucon
        elif self.plot_dict['which_cont'] == 'World':
            self.data = np.nan_to_num(np.true_divide(self.data1, self.pop.reshape(-1, 1)) * self.per)
            self.rdata = self.data1
            self.cntryname = list(k.data_dict['sort_conf']['Country/Region'])

        self.ind, self.pos, self.smoth = smooth_transition(self.data, self.numt, self.topnum)
        self.rind, self.rpos, self.rsmoth = smooth_transition(self.rdata, self.numt, self.topnum)


    def getbardata(self, num):

        self.num = num
        self.numname = self.num % self.intpolnum

        # percapita data
        self.ytop = self.topnum + 1.5
        self.yall = self.smoth.shape[0] - self.smoth[:, self.num]

        self.newtopnum = np.size(np.nonzero(self.yall < self.ytop)[0])
        self.topidx = np.flip(self.ind[:, self.num])[:self.newtopnum]
        self.ynames = [self.cntryname[i] for i in self.topidx]

        self.y = self.smoth.shape[0] - self.smoth[self.topidx, self.num]
        self.xdata = self.data[self.topidx, self.num]

        # raw data
        self.ryall = self.rsmoth.shape[0] - self.rsmoth[:, self.num]

        self.rnewtopnum = np.size(np.nonzero(self.ryall < self.ytop)[0])
        self.rtopidx = np.flip(self.rind[:, self.num])[:self.rnewtopnum]
        self.rynames = [self.cntryname[i] for i in self.rtopidx]

        self.ry = self.rsmoth.shape[0] - self.rsmoth[self.rtopidx, self.num]
        self.rxdata = self.rdata[self.rtopidx, self.num]

        # determine if topidx contains zeros
        testconf = self.data[self.topidx, self.num]
        self.numsigdx = len(np.nonzero(testconf)[0])
        rtestconf = self.rdata[self.rtopidx, self.num]
        self.rnumsigdx = len(np.nonzero(rtestconf)[0])

    def plotbargraph(self, title, rect, rect2, bgimg, adnice):

        # set figure and axis
        fig = plt.figure(figsize=(self.plot_dict['myres_w'] / self.plot_dict['mydpi'],
                                  self.plot_dict['myres_h'] / self.plot_dict['mydpi']),
                         dpi=self.plot_dict['mydpi'], facecolor=self.plot_dict['bkgcolor'])

        self.fig = fig
        if adnice == 1:
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

        coloralpha = 0.4
        self.contcolor = []
        for i, cntry in enumerate(self.ynames):
            # print(cntry)
            cncode = \
            list(self.data_dict['sort_conf'].loc[self.data_dict['sort_conf']['Country/Region'] == cntry, 'code'])[
                0].lower()
            colorname = np.array(
                colors.ColorConverter.to_rgba(list(flagname.loc[flagname['code'] == cncode, 'color'])[0]))
            colorname2 = colorname.copy()
            colorname[-1] = coloralpha
            self.contcolor.append(colorname2)
            self.h[i].set_color(colorname)
            # self.h[i].set_edgecolor(coloredge)
        new_patches = []
        for patch in reversed(self.ax.patches):
            bb = patch.get_bbox()
            self.color = patch.get_facecolor()
            # self.edge = patch.get_edgecolor()
            p_bbox = FancyBboxPatch((bb.xmin, bb.ymin),
                                    abs(bb.width), abs(bb.height),
                                    boxstyle="round, pad=0, rounding_size=0.001",
                                    ec="none", fc=self.color,
                                    mutation_aspect=4,
                                    path_effects=[path_effects.withSimplePatchShadow(offset=(-0.3, -0.3),
                                                                                     shadow_rgbFace=(0, 0, 0),
                                                                                     alpha=0.6)]
                                    )
            patch.remove()
            new_patches.append(p_bbox)
        for patch in new_patches:
            self.ax.add_patch(patch)


        # set axis properties
        if np.mean(self.plot_dict['bkgcolor'][0]) < 0.5:
            self.headcolor = [1,1,1]
            nowcolor = [self.headcolor]*self.numsigdx
        else:
            self.headcolor = [0,0,0]
            nowcolor = self.contcolor[:self.numsigdx]

        self.ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
        self.ax.set_yticks(self.y[:self.numsigdx])
        self.ax.set_yticklabels(self.ynames[:self.numsigdx],
                                fontdict={'fontname': self.plot_dict['myfont'], 'fontweight': 'semibold', 'fontsize': 4.3,
                                          'va': 'center'})
        for ytick, color in zip(self.ax.get_yticklabels(), nowcolor):
            ytick.set_color(color)
        self.ax.tick_params(axis='both', which='both', length=0, pad=2)
        xtickrange = np.round((np.arange(0, self.maxdata, stepgrid))).astype('int')
        self.ax.set_ylim(0.68, self.topnum + 0.5)
        self.ax.set_xticks(xtickrange)
        self.ax.set_xlim(0, 1.2 * xlimit)
        self.ax.set_xticklabels(xtickrange,
                                fontdict={'fontname': self.plot_dict['myfont'], 'fontsize': 4, 'fontweight': 3},
                                color= self.headcolor + [0.6])
        plt.grid(axis='x', c='gray', ls='-', lw=0.4, alpha=0.3)

        #########################################################################################################
        # PLOT GRAPH FOR RAW DATA
        self.rax = self.fig.add_axes(rect2)
        self.rh = self.rax.barh(self.ry, self.rxdata, height=0.7)
        self.rmaxdata = np.max(self.rxdata)

        xlimit = self.rmaxdata * 1.05
        scal = np.array([i * (10 ** (len(str(int(xlimit))) - 2)) for i in [1, 2, 5, 10, 20, 50]])
        ngrid = np.array([xlimit // i for i in scal])
        stepgrid = np.min(scal[ngrid < 10])

        coloralpha = 0.4
        self.rcontcolor = []
        for i, cntry in enumerate(self.rynames):
            # print(cntry)
            cncode = \
                list(self.data_dict['sort_conf'].loc[self.data_dict['sort_conf']['Country/Region'] == cntry, 'code'])[
                    0].lower()
            colorname = np.array(
                colors.ColorConverter.to_rgba(list(flagname.loc[flagname['code'] == cncode, 'color'])[0]))
            # coloredge = np.array(colors.ColorConverter.to_rgba(list(shortname.loc[shortname['Team'] == cntry, 'edge'])[0]))
            colorname2 = colorname.copy()
            colorname[-1] = coloralpha
            self.rcontcolor.append(colorname2)
            self.rh[i].set_color(colorname)
            # self.h[i].set_edgecolor(coloredge)
        new_patches = []
        for patch in reversed(self.rax.patches):
            bb = patch.get_bbox()
            self.color = patch.get_facecolor()
            # self.edge = patch.get_edgecolor()
            p_bbox = FancyBboxPatch((bb.xmin, bb.ymin),
                                    abs(bb.width), abs(bb.height),
                                    boxstyle="round, pad=0, rounding_size=0.001",
                                    ec="none", fc=self.color,
                                    mutation_aspect=4,
                                    path_effects=[path_effects.withSimplePatchShadow(offset=(-0.3, -0.3),
                                                                                     shadow_rgbFace=(0, 0, 0),
                                                                                     alpha=0.6)]
                                    )
            patch.remove()
            new_patches.append(p_bbox)
        for patch in new_patches:
            self.rax.add_patch(patch)

        # set axis properties
        if np.mean(self.plot_dict['bkgcolor'][0]) < 0.5:
            self.headcolor = [1, 1, 1]
            nowcolor = [self.headcolor] * self.rnumsigdx
        else:
            self.headcolor = [0, 0, 0]
            nowcolor = self.rcontcolor[:self.rnumsigdx]

        self.rax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, labelleft=False, labelright=True)
        self.rax.set_yticks(self.ry[:self.rnumsigdx])
        self.rax.set_yticklabels(self.rynames[:self.rnumsigdx],
                                fontdict={'fontname': self.plot_dict['myfont'], 'fontweight': 'semibold', 'fontsize': 4.3,
                                          'va': 'center'})
        for ytick, color in zip(self.rax.get_yticklabels(), nowcolor):
            ytick.set_color(color)
        self.rax.tick_params(axis='both', which='both', length=0, pad=2)
        xtickrange = np.round((np.arange(0, self.rmaxdata, stepgrid))).astype('int')
        self.rax.set_ylim(0.68, self.topnum + 0.5)
        self.rax.set_xticks(xtickrange)
        self.rax.set_xlim(0, 1.2 * xlimit)
        self.rax.set_xticklabels(xtickrange,
                                fontdict={'fontname': self.plot_dict['myfont'], 'fontsize': 4, 'fontweight': 3},
                                color=self.headcolor + [0.6])
        plt.grid(axis='x', c='gray', ls='-', lw=0.4, alpha=0.3)


    def add_items(self, zoomscalar = 0.018):

        # insert flag and annotate the number of cases at the right of it
        self.zoomscalar = zoomscalar
        imagedist = self.maxdata / 20
        for i, cntry in enumerate(self.ynames):
            if i < self.numsigdx:
                cncode = list(self.data_dict['sort_conf'].loc[self.data_dict['sort_conf']['Country/Region'] == cntry, 'code'])[0].lower()
                c = f'worldflag/{cncode}.png'
                offset_image(self.xdata[i] + imagedist, self.y[i], c, self.ax, self.fig, self.zoomscalar)

        for py, px in enumerate(self.xdata):
            text1 = np.round(px, 1)
            contname = list(self.data_dict['sort_conf'].loc[
                                self.data_dict['sort_conf']['Country/Region'] == self.ynames[py], 'Continent_Name'])[
                0].upper()
            if contname == 'NORTH AMERICA':
                contname = 'N. AMERICA'
            elif contname == 'SOUTH AMERICA':
                contname = 'S. AMERICA'
            if py < self.numsigdx:
                self.ax.annotate(text1, xy=(px + imagedist + self.maxdata / 18, self.y[py]),
                                 fontname=self.plot_dict['myfont'],
                                 fontsize=4.2, fontweight='semibold', fontstyle='normal', va='center', ha='left',
                                 color=self.headcolor)
                if self.plot_dict['which_cont'] == 'World':
                    if (len(contname) > 9):
                        if ((px / self.maxdata) > 0.18):
                            self.ax.annotate(contname, xy=(px - self.maxdata / 130, self.y[py]), color='w',
                                             fontweight='semibold',
                                             fontname=self.plot_dict['myfont'], fontsize=4.2, fontstyle='normal',
                                             va='center', ha='right')
                    else:
                        if ((px / self.maxdata) > 0.12):
                            self.ax.annotate(contname, xy=(px - self.maxdata / 130, self.y[py]), color='w',
                                             fontweight='semibold',
                                             fontname=self.plot_dict['myfont'], fontsize=4.2, fontstyle='normal',
                                             va='center', ha='right')

        #######################################################################################
        ## FOR  RAW DATA
        imagedist = self.rmaxdata / 20
        for i, cntry in enumerate(self.rynames):
            if i < self.rnumsigdx:
                cncode = \
                list(self.data_dict['sort_conf'].loc[self.data_dict['sort_conf']['Country/Region'] == cntry, 'code'])[
                    0].lower()
                c = f'worldflag/{cncode}.png'
                offset_image(self.rxdata[i] + imagedist, self.ry[i], c, self.rax, self.fig, self.zoomscalar)

        for py, px in enumerate(self.rxdata):
            text1 = int2str(str(int(px)))
            contname = list(self.data_dict['sort_conf'].loc[
                                self.data_dict['sort_conf']['Country/Region'] == self.rynames[py], 'Continent_Name'])[
                0].upper()
            if contname == 'NORTH AMERICA':
                contname = 'N. AMERICA'
            elif contname == 'SOUTH AMERICA':
                contname = 'S. AMERICA'
            if py < self.rnumsigdx:
                self.rax.annotate(text1, xy=(px + imagedist + self.rmaxdata / 18, self.ry[py]),
                                 fontname=self.plot_dict['myfont'],
                                 fontsize=4.2, fontweight='semibold', fontstyle='normal', va='center', ha='right',
                                 color=self.headcolor)
                if self.plot_dict['which_cont'] == 'World':
                    if (len(contname) > 9):
                        if ((px / self.rmaxdata) > 0.18):
                            self.rax.annotate(contname, xy=(px - self.rmaxdata / 130, self.ry[py]), color='w',
                                             fontweight='semibold',
                                             fontname=self.plot_dict['myfont'], fontsize=4.2, fontstyle='normal',
                                             va='center', ha='left')
                    else:
                        if ((px / self.rmaxdata) > 0.12):
                            self.rax.annotate(contname, xy=(px - self.rmaxdata / 130, self.ry[py]), color='w',
                                             fontweight='semibold',
                                             fontname=self.plot_dict['myfont'], fontsize=4.2, fontstyle='normal',
                                             va='center', ha='left')


    def add_text(self, sidetop, rect1=[0.44, 0.9, 0.12, 0.1]):
        self.ax1 = self.fig.add_axes(rect1)

        self.ax1.text(0.5, 0.5, '{} {}'.format(self.day[self.num], self.month[self.num]), transform=self.ax1.transAxes,fontweight='semibold',
                     verticalalignment='center', horizontalalignment='center', color=self.headcolor, fontsize=18, fontname=self.plot_dict['myfont'])

        tworld = self.wall[self.num]
        xs, ys, ystp = 0.7, 0.2, 0.08
        bigfont, smallfont = 12, 10

        self.ax.text(xs, ys, 'WORLD', transform=self.ax.transAxes,
                      va='center', ha='center', color=self.headcolor, fontweight='semibold', fontsize=bigfont,
                      fontname=self.plot_dict['myfont'])
        self.ax.text(xs, ys - ystp, '{:2.1f}\nper 1M pop'.format(tworld), transform=self.ax.transAxes,
                      va='center', ha='center',
                      color=self.headcolor, fontsize=smallfont, fontweight='semibold',
                      fontname=self.plot_dict['myfont'])

        tworldsum = np.sum(self.data_w[:, self.num])
        text = str(int(tworldsum))
        text1 = int2str(text)
        self.rax.text(1-xs, ys, 'WORLD', transform=self.rax.transAxes,
                      va='center', ha='center', color=self.headcolor, fontweight='semibold', fontsize=bigfont,
                      fontname=self.plot_dict['myfont'])
        self.rax.text(1-xs, ys - ystp, '{}'.format(text1), transform=self.rax.transAxes, va='center', ha='center',
                      color=self.headcolor, fontsize=smallfont, fontweight='semibold',
                      fontname=self.plot_dict['myfont'])

        self.ax1.set_axis_off()

    def set_finalaxis(self):
        self.ax.invert_yaxis()
        self.ax.set_frame_on(False)
        self.rax.invert_yaxis()
        self.rax.invert_xaxis()
        self.rax.set_frame_on(False)

    def saveplot(self):
        # write to video
        figname = 'D:/YoutubeChannel/CoronaVirus/CoronaImages/Bar2{}/{}{}-{}.png'.format(self.plot_dict['which_cont'], self.month[self.num],
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
    rect2 = [0.01, 0.05, 0.38, 0.85]
    rect = [0.61, 0.05, 0.38, 0.85]
    adnice = 0
    k.plotbargraph(title, rect, rect2, bgimg, adnice)
    zoomscalar = 0.018
    k.add_items(zoomscalar)
    sidetop = 3
    k.add_text(sidetop)
    k.set_finalaxis()
    k.saveplot()

    return f'done with {num}'

################################################################################################################
where_to_start = 0
startnow = 0
data_dict = dict(sort_conf=sort_conf, sort_death=sort_death, sort_rec=sort_rec, conf=conf, deat=deat, recd=recd)
myfont = 'Consolas'
bkgcolor = np.array([1,1,1,1])-0.9
bkgcolor[-1] = 0.95
which_cont = 'Africa'
plot_dict = dict(where_to_start=where_to_start, myfont=myfont, bkgcolor=bkgcolor, myres_w=1920,
                 myres_h=1080, mydpi=300, which_cont=which_cont)

datatype='Total'
k = Create_Animation(data_dict, datatype)
Create_Animation.set_variables(plot_dict)
intpolnum = 90
topnum = 20
numt = 60
k.intpoldata(intpolnum, topnum, numt)
totalnum = intpolnum * (len(sort_conf.columns[6:]) - where_to_start)
startnum = 0  #intpolnum * startnow

num = totalnum-1
plotgraph(num)

#################
#
# if __name__ == '__main__':
#
#      with concurrent.futures.ProcessPoolExecutor() as executor:
#         num = range(startnum, totalnum, 400)
#         results = executor.map(plotgraph, num)
#
#         for result in results:
#             print(result)
#



