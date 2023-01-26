## CREATING ANIMATION FOR MAP AFRICA

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as color
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import cv2 as cv
import geopandas as gpd
import matplotlib.pyplot as plt
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


def offset_image(xcoord, ycoord, name, ax, zoomscalar, xres=600, yres=400):
    img1 = get_image(name)
    img = cv.resize(img1, (xres, yres), interpolation=cv.INTER_AREA)
    im = OffsetImage(img, zoom=zoomscalar, transform=ax.transAxes)
    im.image.axes = ax
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


# Maps values to a bin.
# The mapped values must start at 0 and end at 1.
def bin_mapping(x):
    # mybin = [100,200,1000,5000,10000,20000,50000,100000,500000, 1000000, max_cases]
    for idx, bound in enumerate(mybin):
        if x < bound:
            return idx / (len(mybin) - 1.0)


def my_annotate(x, y, x1, y1, label):
    plt.annotate(label, xy=(x, y), xytext=(x + x1, y + y1),
                 arrowprops=dict(facecolor='gray', arrowstyle='-'),
                 fontsize=7, verticalalignment='center', horizontalalignment='center', color='black', fontname=myfont)


def changestate(confirmed, death, rec, chdict):
    for t in np.arange(len(chdict)):
        confirmed.loc[confirmed['Country/Region'] == chdict[t][0], 'Country/Region'] = chdict[t][1]
        death.loc[death['Country/Region'] == chdict[t][0], 'Country/Region'] = chdict[t][1]
        rec.loc[rec['Country/Region'] == chdict[t][0], 'Country/Region'] = chdict[t][1]
    return confirmed, death, rec


# ## Create Parent Class
# #### The parent class will contain class variables comprising of all cleaned data, data accessorys and plot variables

class Create_Animation():
    plot_dict = None
    acc_dict = None


    def __init__(self, data_dict):
        self.data_dict = data_dict

    @classmethod
    def set_variables(cls, plot_dict):
        cls.plot_dict = plot_dict

    def intpoldata(self, intpolnum):
        self.intpolnum = intpolnum
        self.days = self.data_dict['sort_conf'].columns[4:]
        self.conf_all = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['sort_conf'],
                                   self.intpolnum)
        self.death_all = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['sort_death'],
                                    self.intpolnum)
        self.rec_all = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['sort_rec'],
                                  self.intpolnum)
        self.conf_w = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['conf'], self.intpolnum)
        self.death_w = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['deat'], self.intpolnum)
        self.rec_w = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['recd'], self.intpolnum)
        self.day, self.month = daylist(self.days, self.plot_dict['where_to_start'], self.intpolnum)

    def getmapdata(self, num, datatype):
        self.num = num
        self.numname = self.num % self.intpolnum
        self.act_all = self.conf_all - (self.death_all + self.rec_all)

        if datatype == 'Active':
            nowdata = self.act_all[:, self.num]
        elif datatype == 'Total':
            nowdata = self.conf_all[:, self.num]
        elif datatype == 'Deaths':
            nowdata = self.death_all[:, self.num]
        elif datatype == 'Recovered':
            nowdata = self.rec_all[:, self.num]


        self.nodf = self.data_dict['sort_conf'][['code', 'Country/Region']].copy()
        self.nodf['data'] = nowdata
        self.nodf = self.nodf.rename(columns={'code': 'ISO2'})
        self.nowdf = pd.merge(self.plot_dict['world_copy'].copy(), self.nodf, on='ISO2', how='inner')
        self.nowvalues = self.nowdf['data'].values
        self.nowdf['databin'] = self.nowdf['data'].apply(bin_mapping)

        conf = self.data_dict['sort_conf'][['code', 'Country/Region']].copy()
        conf['data'] = self.conf_all[:, self.num]
        conf = conf.rename(columns={'code': 'ISO2'})
        self.confdf = pd.merge(self.plot_dict['world_copy'].copy(), conf, on='ISO2', how='inner')

        act = self.data_dict['sort_conf'][['code', 'Country/Region']].copy()
        act['data'] = self.act_all[:, self.num]
        act = act.rename(columns={'code': 'ISO2'})
        self.actdf = pd.merge(self.plot_dict['world_copy'].copy(), act, on='ISO2', how='inner')

        dedf = self.data_dict['sort_death'][['code', 'Lat', 'Long']].copy()
        dedf['data'] = self.death_all[:, self.num]
        dedf = dedf.rename(columns={'code': 'ISO2'})
        self.deathdf = pd.merge(self.plot_dict['world_copy'].copy(), dedf, on='ISO2', how='inner')

        recf = self.data_dict['sort_rec'][['code', 'Lat', 'Long']].copy()
        recf['data'] = self.rec_all[:, self.num]
        recf = recf.rename(columns={'code': 'ISO2'})
        self.recdf = pd.merge(self.plot_dict['world_copy'].copy(), recf, on='ISO2', how='inner')

    def plotmap(self, title, Legend_elements, cmap, world_points, rect=[0, 0.05, 1, 1]):
        # open figure
        fig = plt.figure(figsize=(self.plot_dict['myres_w'] / self.plot_dict['mydpi'],
                                  self.plot_dict['myres_h'] / self.plot_dict['mydpi']),
                         dpi=self.plot_dict['mydpi'], facecolor=self.plot_dict['bkgcolor'])
        self.fig = fig
        self.ax = self.fig.add_axes(rect)

        self.nowdf.plot(column=self.nowdf['databin'], ax=self.ax, cmap=cmap,
                        legend=False, edgecolor='gray', linewidth=0.2, vmin=0, vmax=1)

        self.ax.set_frame_on(False)
        self.ax.set_axis_off()
        # self.ax.set_title(title, color='red',
        #                   fontdict={'fontsize': 15, 'fontname': self.plot_dict['myfont'], 'ha': 'center',
        #                             'fontweight': 'semibold', 'fontstretch': 'semi-condensed'})
        self.ax.legend(handles=Legend_elements, bbox_to_anchor=(0.2, 0.6), fontsize=6,
                       markerscale=0.7, frameon=False)

        # add country names
        texts = []
        wgt = 'semibold'
        for x, y, iso, iso3, label in zip(world_points.geometry.x, world_points.geometry.y, world_points["ISO2"], world_points["ISO3"],
                                    world_points["NAME"]):
            textcolor = 'k'
            self.iso = iso
            if (self.nowdf.loc[self.nowdf['ISO2'] == iso, 'databin'] > (6/11)).all():
                textcolor = 'w'
            else:
                textcolor = 'k'
            if iso == 'US':
                x, y = x + 10, y - 6
            elif iso == 'CA':
                x, y = x - 10, y - 3
            elif iso3 == 'NOR':
                x,y = x-5, y-5

            if world_points.loc[world_points['ISO2'] == iso, 'area'].values[0] > 278.5:
                texts.append(plt.text(x, y, label, fontsize=7, verticalalignment='center', fontweight=wgt, fontstretch='semi-condensed',
                                      horizontalalignment='center', color=textcolor, fontname=self.plot_dict['myfont']))
            elif (world_points.loc[world_points['ISO2'] == iso, 'area'].values[0] > 100) and (world_points.loc[world_points['ISO2'] == iso, 'area'].values[0] < 278.5):
                texts.append(plt.text(x, y, iso3, fontsize=5, verticalalignment='center', fontweight=wgt,fontstretch='semi-condensed',
                                      horizontalalignment='center', color=textcolor, fontname=self.plot_dict['myfont']))
            elif (world_points.loc[world_points['ISO2'] == iso, 'area'].values[0] > 30) and (world_points.loc[world_points['ISO2'] == iso, 'area'].values[0] < 100):
                texts.append(plt.text(x, y, iso3, fontsize=4, verticalalignment='center', fontweight=wgt,fontstretch='semi-condensed',
                                      horizontalalignment='center', color=textcolor, fontname=self.plot_dict['myfont']))
            elif (world_points.loc[world_points['ISO2'] == iso, 'area'].values[0] > 10) and (world_points.loc[world_points['ISO2'] == iso, 'area'].values[0] < 30):
                if iso == 'SJ':
                    iso = ''
                texts.append(plt.text(x, y, iso, fontsize=4, verticalalignment='center', fontweight=wgt,fontstretch='semi-condensed',
                                      horizontalalignment='center', color=textcolor, fontname=self.plot_dict['myfont']))



    def add_basictext(self, topnum, rect2=[0.35, -0.03, 0.3, 0.32]):

        zoomscalar = 0.03
        wgt = 'semibold'
        # get the top countries

        topidx = np.flip(np.argsort(self.nowvalues))[:topnum]

        # determine if topidx contains zeros
        testconf = self.nowvalues[topidx]
        numsigdx = len(np.nonzero(testconf)[0])
        topidx = topidx[range(numsigdx)]

        ynames = list(self.nowdf.loc[topidx, 'NAME'])
        xconf = self.confdf['data'].values[topidx]
        xact = self.actdf['data'].values[topidx]
        xdeath = self.deathdf['data'].values[topidx]
        xrec = self.recdf['data'].values[topidx]

        self.ax.text(0.89, 1, '{} {}'.format(self.day[self.num], self.month[self.num]), transform=self.ax.transAxes,fontweight=wgt,
                     verticalalignment='bottom', horizontalalignment='center', color='black', fontsize=20,
                     fontname=self.plot_dict['myfont'])
        self.ax.text(0.58, 0.13, 'Top {} Countries'.format(topnum), transform=self.ax.transAxes,fontweight=wgt,
                     verticalalignment='bottom', horizontalalignment='center', color='red', fontsize=10,
                     fontname=self.plot_dict['myfont'])

        #     # add flags
        nowfont = 6
        nm, fl, cs, at, dth, rc = 0.2, 0.35, 0.55, 0.8, 1.05, 1.25
        ax2 = self.fig.add_axes(rect2, anchor='SE', zorder=-1)
        ax2.text((fl + nm) / 2, 0.93, 'States', transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(cs, 0.93, 'Total(#)', transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(at, 0.93, 'Active(#)', transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(dth, 0.93, 'Dths(%)', transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(rc, 0.93, 'Rcvd(%)', transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=wgt,
                 fontname=self.plot_dict['myfont'])

        ystart = 0.8
        ystp = 0.2
        for i, cntry in enumerate(ynames):
            cncode = list(flagname.loc[flagname['Country'] == cntry, 'code'])[0]
            c = 'worldflag/' + cncode + '.png'
            offset_image(fl, ystart - i * ystp, c, ax2, zoomscalar)
            label = list(self.plot_dict['world_copy'].loc[self.plot_dict['world_copy']['ISO2'] == cncode.upper(), 'ISO3'])[0]
            ax2.text(nm, ystart - i * ystp, '{}'.format(label), transform=ax2.transAxes,
                     verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=wgt,
                     fontname=self.plot_dict['myfont'])
            ax2.text(cs, ystart - i * ystp, '{:8d}'.format(int(xconf[i])), transform=ax2.transAxes,
                     verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=wgt,
                     fontname=self.plot_dict['myfont'])
            ax2.text(at, ystart - i * ystp, '{:8d}'.format(int(xact[i])), transform=ax2.transAxes,
                     verticalalignment='center', horizontalalignment='center', color='b', fontsize=nowfont,fontweight=wgt,
                     fontname=self.plot_dict['myfont'])
            ax2.text(dth, ystart - i * ystp, '{:2.2f}'.format((xdeath[i] / xconf[i]) * 100),
                     transform=ax2.transAxes,
                     verticalalignment='center', horizontalalignment='center', color='red', fontsize=nowfont,fontweight=wgt,
                     fontname=self.plot_dict['myfont'])
            ax2.text(rc, ystart - i * ystp, '{:2.2f}'.format((xrec[i] / xconf[i]) * 100), transform=ax2.transAxes,
                     verticalalignment='center', horizontalalignment='center', color='g', fontsize=nowfont,fontweight=wgt,
                     fontname=self.plot_dict['myfont'])

        ax2.text((fl + nm) / 2, ystart - topnum * ystp, '{}'.format('WORLD'), transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='black', fontsize=nowfont + 2,fontweight=wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(cs, ystart - topnum * ystp, '{:8d}'.format(int(np.sum(self.conf_all[:, self.num]))),
                 transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='black', fontsize=nowfont,fontweight=wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(at, ystart - topnum * ystp, '{:8d}'.format(int(np.sum(self.act_all[:, self.num]))),
                 transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='b', fontsize=nowfont,fontweight=wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(dth, ystart - topnum * ystp,
                 '{:2.2f}'.format((np.sum(self.death_all[:, self.num]) / np.sum(self.conf_all[:, self.num]))*100),
                 transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='r', fontsize=nowfont,fontweight=wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(rc, ystart - topnum * ystp,
                 '{:2.2f}'.format((np.sum(self.rec_all[:, self.num]) / np.sum(self.conf_all[:, self.num]))*100),
                 transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='g', fontsize=nowfont,fontweight=wgt,
                 fontname=self.plot_dict['myfont'])

        ax2.set_axis_off()

    def add_image(self, corimg, rect1=[0.05, 0.85, 0.1, 0.1]):
        # inset coronavirus image
        ax1 = self.fig.add_axes(rect1, anchor='NW', zorder=-1)
        offset_image2(0.5, 0.5, corimg, ax1, 0.07, 500, 500)
        ax1.set_axis_off()

    # #     # add situational reports box
    def add_report(self, report):

        xst = 0.55
        yst = -0.1
        wgt = 'semibold'
        edcolor = 'black'
        talpha = 0
        textlim = 50
        tcolor = 'red'
        dayidx = int(self.num // self.intpolnum) + self.plot_dict['where_to_start']
        rpnum = int(type(report.iloc[dayidx, 1]) == str)
        if rpnum == 1:
            if self.numname > int(self.intpolnum * 0.1) and self.numname < int(self.intpolnum * 0.9):
                text1 = normalize('NFKD', report.iloc[dayidx, 1])
                tx = text1.split()
                cumtx = np.cumsum(np.array([len(i) + 1 for i in tx]))
                self.newtxt = ''
                n = 1

                while len(text1) > len(self.newtxt):
                    lastindx = np.int(np.nonzero(cumtx < textlim * n)[0][-1])
                    if n == 1:
                        self.newtxt = " ".join(tx[: lastindx + 1])
                    else:
                        self.newtxt = self.newtxt + '\n' + " ".join(tx[previndx + 1: lastindx + 1])

                    previndx = lastindx
                    n += 1

                self.ax.text(xst, yst, '{}'.format(self.newtxt), transform=self.ax.transAxes,
                             verticalalignment='center', horizontalalignment='right', color=tcolor, fontsize=7,fontweight=wgt,
                             fontname=self.plot_dict['myfont'],
                             bbox=dict(fc='w', ec=edcolor, pad=5, alpha=talpha))

    def saveplot(self):
        # write to video
        figname = 'D:/YoutubeChannel/CoronaVirus/CoronaImages/World/{}{}-{}.png'.format(self.month[self.num],
                                                                                         self.day[self.num], self.numname)
        self.fig.savefig(figname, dpi=self.fig.get_dpi(), facecolor=self.fig.get_facecolor())
        plt.cla()
        self.fig.clf()
        plt.close('all')
        gc.collect()


def plotgraph(num):
    k.getmapdata(num, datatype)

    title = 'Global Active Cases of COVID-19'
    topnum = 3

    k.plotmap(title, Legend_elements, cmap, world_points)
    k.add_basictext(topnum)
    # k.add_image(corimg)
    # k.add_report(report)
    k.saveplot()

    return f'done with {num}'


################################################################################################################
# ## Import and cleaning data

compath = "C:/Users/copok/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
dpath = "C:/Users/copok/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
rpath = "C:/Users/copok/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"
confirmed = pd.read_csv(compath)
death = pd.read_csv(dpath)
rec = pd.read_csv(rpath)
print('The shapes of confirmed, death and recoverd are {}, {} and {}'.format(confirmed.shape, death.shape, rec.shape))

confirmed.loc[confirmed['Province/State'] == 'Greenland', 'Country/Region'] = 'Greenland'
death.loc[death['Province/State'] == 'Greenland', 'Country/Region'] = 'Greenland'
rec.loc[rec['Province/State'] == 'Greenland', 'Country/Region'] = 'Greenland'

chdict = [['Korea, South','South Korea'], ['Taiwan*','Taiwan'], ['Burma', 'Myanmar']]
confirmed, death, rec = changestate(confirmed, death, rec, chdict)

conf = cleandata(confirmed)
deat = cleandata(death)
recd = cleandata(rec)
print('The shapes of confirmed, death and recovered are {}, {} and {}'.format(conf.shape, deat.shape, recd.shape))

# download us codes and choose the continent
uscode = pd.read_csv("D:/YoutubeChannel/CoronaVirus/worldcodes2.csv")
# uscode = uscode[uscode.Continent_Name == 'Africa']
uscode = uscode.rename(columns={'Country': 'Country/Region'})
uscode.loc[uscode['Country/Region'] == 'Namibia', 'code'] = 'NA'
sort_conf = pd.merge(uscode[['code', 'Country/Region']], conf, on='Country/Region', how='inner')
sort_death = pd.merge(uscode[['code', 'Country/Region']], deat, on='Country/Region', how='inner')
sort_rec = pd.merge(uscode[['code', 'Country/Region']], recd, on='Country/Region', how='inner')
print('The shapes of confirmed, death and recovered are {}, {} and {}'.format(sort_conf.shape,
                                                                              sort_death.shape, sort_rec.shape))

namematch = [['KP', 'North Korea']]
mainlist = []
for ml in namematch:
    r = ml
    for i in np.arange(sort_conf.shape[1] - 2):
        r.append(0)
    mainlist.append(r)

ld = pd.DataFrame(mainlist, columns=list(sort_conf.columns))
sort_conf = sort_conf.append(ld)
sort_death = sort_death.append(ld)
sort_rec = sort_rec.append(ld)

# ## Import Accessory Data such as Flags, Images and world map
# download situation reports and flags
report = pd.read_excel("D:/YoutubeChannel/CoronaVirus/SituationReports.xlsx")
flagpath = "D:/YoutubeChannel/CoronaVirus/worldflags.csv"
flagname = pd.read_csv(flagpath)

# get coronavirus image
corimg = get_image('coronavirus.png')

shapepath = "D:/YoutubeChannel/CoronaVirus/worldshapefile/TM_WORLD_BORDERS-0.3.shp"
# shapepath = "D:/YoutubeChannel/CoronaVirus/worldshapefile/ne_50m_admin_0_countries.shp"
# world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
world = gpd.read_file(shapepath)
world = world[world.NAME != "Antarctica"]
'Iran (Islamic Republic of)'

world_copy = world.copy()
world_copy.loc[world_copy['NAME'] == 'Korea, Republic of', 'NAME'] = 'South Korea'
world_copy.loc[world_copy['NAME'] == 'Iran (Islamic Republic of)', 'NAME'] = 'Iran'

world_copy["center"] = world_copy["geometry"].centroid
world_points = world_copy.copy()
world_points = world_points.set_geometry("center")
world_points['area'] = world_copy['geometry'].area

# ## Setting Plot Variables
# SET LEGEND
datatype = 'Recovered'
if datatype == 'Active' or datatype == 'Total':
    conlist = list(sort_conf.iloc[:, 3:].max())
    max_cases = max(conlist) + 100
    print('The maximum cases is {}'.format(max_cases))
    mybin = [100,200,1000,5000,10000,20000,50000,100000,500000, 1000000, max_cases]
elif datatype == 'Recovered':
    conlist = list(sort_rec.iloc[:, 3:].max())
    max_cases = max(conlist) + 100
    print('The maximum cases is {}'.format(max_cases))
    mybin = [100, 200, 1000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, max_cases]
elif datatype == 'Deaths':
    conlist = list(sort_death.iloc[:, 3:].max())
    max_cases = max(conlist) + 100
    print('The maximum cases is {}'.format(max_cases))
    mybin = [10, 100, 200, 500, 1000, 5000, 10000, 20000, 50000, 100000, max_cases]

lglist = []
for i in np.arange(len(mybin)):
    if i == 0:
        lglist.append('<{}'.format(mybin[i]))
    elif i == len(mybin) - 1:
        lglist.append('{}>'.format(mybin[i - 1] + 1))
    else:
        lglist.append('{}-{}'.format(mybin[i - 1] + 1, mybin[i]))

# SET COLOR MAP
c1 = mpl.cm.get_cmap('RdYlBu_r', 20)
c2 = mpl.cm.get_cmap('gist_heat_r', 20)
colorlist = np.vstack((c1(range(4, 20, 2)), c2(range(15, 20, 2))))
colorlist[-2] = colorlist[-1]
colorlist[-1] = [0.2, 0, 0.4, 1]
bin_labels = [idx / (len(mybin) - 1.0) for idx in range(len(mybin))]
# Create the custom color map
cmap = color.LinearSegmentedColormap.from_list('mycmap',
                                               [(lbl, color) for lbl, color in zip(bin_labels, colorlist)])

mksize = 10
bkgcolor = [1,1,1,0.8]
myfont = 'Consolas'
Legend_elements = []
for i in np.arange(len(mybin)):
    Legend_elements.append(Line2D([0], [0], marker='o', markersize=mksize, color=bkgcolor,
                                  markerfacecolor=colorlist[i], label=lglist[i]))

where_to_start = 134
data_dict = dict(sort_conf=sort_conf, sort_death=sort_death, sort_rec=sort_rec, conf=conf, deat=deat, recd=recd)
# myfont = 'Consolas'
# bkgcolor = 'lightgray'
plot_dict = dict(where_to_start=where_to_start, world_copy=world_copy, myfont=myfont, bkgcolor=bkgcolor, myres_w=1920,
                 myres_h=1080, mydpi=300)

k = Create_Animation(data_dict)
Create_Animation.set_variables(plot_dict)
intpolnum = 90
k.intpoldata(intpolnum)
totalnum = intpolnum * (len(sort_conf.columns[4:]) - where_to_start)
num = 0
plotgraph(num)


# if __name__ == '__main__':
#     num = range(0, totalnum)
#
#     with concurrent.futures.ProcessPoolExecutor() as executor:
#         results = executor.map(plotgraph, num)
#
#         for result in results:
#             print(result)
