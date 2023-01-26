## CREATING TIME-BASED ANIMATION TO ILLUSTRATE SPREAD OF COVID-19 FOR ACROSS THE UNITED STATES USING DATA FROM WORLDOMETER

import numpy as np
import pandas as pd
import matplotlib.colors as color
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import cv2 as cv
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import gc
print('Libraries Imported')


#
def cleandata(data):
    p1 = data.groupby('Province_State').first()
    p1 = p1[['Lat', 'Long_']]
    p2 = data.groupby('Province_State').sum()
    return p1.join(p2.iloc[:, 5:]).reset_index()


def get_image(name):
    path = "G:/YoutubeChannel/CoronaVirus/{}".format(name)
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
    # mybin = bin = [10,20,50,100,500,1000,5000,10000,50000,100000, max_cases]
    for idx, bound in enumerate(mybin):
        if x < bound:
            return idx / (len(mybin) - 1.0)


def my_annotate(x, y, x1, y1, label):
    plt.annotate(label, xy=(x, y), xytext=(x + x1, y + y1),
                 arrowprops=dict(facecolor='gray', arrowstyle='-'),
                 fontsize=7, verticalalignment='center', horizontalalignment='center', color='black', fontname=myfont)


# ## Import and cleaning data
compath = "C:/Users/copok/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"
dpath = "C:/Users/copok/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv"
# rpath = "D:/YoutubeChannel/CoronaVirus/time_series_covid19_recovered_global.csv"
confirmed = pd.read_csv(compath)
death = pd.read_csv(dpath)
# rec = pd.read_csv(rpath)
confirmed = confirmed[(confirmed['iso2']== 'US') & (confirmed['Lat'] > 0)]
death = death[(death['iso2']== 'US') & (death['Lat'] > 0)]
print('The shapes of confirmed, death are {} and {}'.format(confirmed.shape, death.shape))


conf = cleandata(confirmed)
deat = cleandata(death)
print('The shapes of confirmed, death are {} and {}'.format(conf.shape, deat.shape))

# download us codes and choose the continent
uscode = pd.read_csv("G:/YoutubeChannel/CoronaVirus/codes.csv")
uscode = uscode.rename(columns={'state': 'Province_State'})
sort_conf = pd.merge(uscode[['code', 'Province_State']], conf, on='Province_State', how='inner')
sort_death = pd.merge(uscode[['code', 'Province_State']], deat, on='Province_State', how='inner')
# sort_rec = pd.merge(uscode[['code', 'Province_State']], recd, on='Province_State', how='inner')
print('The shapes of confirmed, death are {} and {}'.format(sort_conf.shape, sort_death.shape))
# conlist = list(sort_conf.iloc[:, 3:].max())
# max_cases = max(conlist) + 100
# print('The maximum cases is {}'.format(max_cases))

# ## Import Accessory Data such as Flags, Images and world map
# download situation reports and flags
report = pd.read_excel("G:/YoutubeChannel/CoronaVirus/SituationReports.xlsx")
flagpath = "G:/YoutubeChannel/CoronaVirus/statesflags.csv"
flagname = pd.read_csv(flagpath)

# get coronavirus image
corimg = get_image('coronavirus.png')

shapepath = "G:/YoutubeChannel/CoronaVirus/USshapefile/cb_2018_us_state_500k.shp"
usa = gpd.read_file(shapepath)
usa_flt = usa.loc[usa['GEOID'].astype('float')<59, :]

usa_copy = usa_flt.copy()
alka = usa_flt.loc[usa['NAME']=='Alaska', 'geometry'].scale(xfact=0.4, yfact=0.4, zfact=1.0, origin='center')
alka = alka.translate(xoff= -55, yoff= -35, zoff=0.0)
usa_copy.loc[usa['NAME']=='Alaska', 'geometry'] = alka

haw = usa_flt.loc[usa['NAME']=='Hawaii', 'geometry'].translate(xoff= 55, yoff= 3, zoff=0.0)
usa_copy.loc[usa['NAME']=='Hawaii', 'geometry'] = haw

usa_copy["center"] = usa_copy["geometry"].centroid
usa_points = usa_copy.copy()
usa_points = usa_points.set_geometry("center")
usa_points['area'] = usa_copy['geometry'].area

# ## Setting Plot Variables
# SET LEGEND
datatype = 'Deaths'
if datatype == 'Total':
    conlist = list(sort_conf.iloc[:, 4:].max())
    max_cases = max(conlist) + 100
    print('The maximum cases is {}'.format(max_cases))
    mybin = [100, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 300000, max_cases]
elif datatype == 'Deaths':
    conlist = list(sort_death.iloc[:, 5:].max())
    max_cases = max(conlist) + 100
    print('The maximum cases is {}'.format(max_cases))
    mybin = [10, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 30000, max_cases]


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
bkgcolor = 'snow'
myfont = 'Consolas'
Legend_elements = []
for i in np.arange(len(mybin)):
    Legend_elements.append(Line2D([0], [0], marker='o', markersize=mksize, color=bkgcolor,
                                  markerfacecolor=colorlist[i], label=lglist[i]))


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
        # self.rec_all = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['sort_rec'],
        #                           self.intpolnum)
        self.conf_w = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['conf'], self.intpolnum)
        self.death_w = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['deat'], self.intpolnum)
        # self.rec_w = all_intpol(self.days, self.plot_dict['where_to_start'], self.data_dict['recd'], self.intpolnum)
        self.day, self.month = daylist(self.days, self.plot_dict['where_to_start'], self.intpolnum)

    def getmapdata(self, num, datatype):
        self.num = num
        self.numname = self.num % self.intpolnum

        if datatype == 'Total':
            nowdata = self.conf_all[:, self.num]
        elif datatype == 'Deaths':
            nowdata = self.death_all[:, self.num]

        self.nodf = self.data_dict['sort_conf'][['code', 'Province_State']].copy()
        self.nodf['data'] = nowdata
        self.nodf = self.nodf.rename(columns={'code': 'STUSPS'})
        self.nowdf = pd.merge(self.plot_dict['usa_copy'].copy(), self.nodf, on='STUSPS', how='inner')
        self.nowvalues = self.nowdf['data'].values
        self.nowdf['databin'] = self.nowdf['data'].apply(bin_mapping)

        conf = self.data_dict['sort_death'][['code', 'Lat', 'Long_']].copy()
        conf['data'] = self.conf_all[:, self.num]
        conf = conf.rename(columns={'code': 'STUSPS'})
        self.confdf = pd.merge(self.plot_dict['usa_copy'].copy(), conf, on='STUSPS', how='inner')

        dedf = self.data_dict['sort_death'][['code', 'Lat', 'Long_']].copy()
        dedf['data'] = self.death_all[:, self.num]
        dedf = dedf.rename(columns={'code': 'STUSPS'})
        self.deathdf = pd.merge(self.plot_dict['usa_copy'].copy(), dedf, on='STUSPS', how='inner')


    def plotmap(self, title, Legend_elements, cmap, usa_points, rect=[0.1, 0, 0.85, 1]):
        # open figure
        fig = plt.figure(figsize=(self.plot_dict['myres_w'] / self.plot_dict['mydpi'],
                                  self.plot_dict['myres_h'] / self.plot_dict['mydpi']),
                         dpi=self.plot_dict['mydpi'], facecolor=self.plot_dict['bkgcolor'])
        self.fig = fig
        self.ax = self.fig.add_axes(rect)
        self.wgt ='bold'

        self.nowdf.plot(column=self.nowdf['databin'], ax=self.ax, cmap=cmap,
                        legend=False, edgecolor='gray', linewidth=0.5, vmin=0, vmax=1)

        self.ax.set_xlim(-130, -60)
        self.ax.set_frame_on(False)
        self.ax.set_axis_off()
        # self.ax.set_title(title, color='red',
        #                   fontdict={'fontsize': 25, 'fontname': self.plot_dict['myfont'], 'fontweight': 3,
        #                             'ha': 'center'})
        L = self.ax.legend(handles=Legend_elements, bbox_to_anchor=(0.1, 0.6), fontsize=6,
                       markerscale=0.7, frameon=False)
        plt.setp(L.texts, family=self.plot_dict['myfont'])

        # add country names
        texts = []
        smstate = dict(NJ=[3, -1], MA=[4, 1], DC=[3, -4], DE=[3, -2], CT=[4, -2], RI=[4, 0], HI=[0, 0])
        for x, y, label in zip(usa_points.geometry.x, usa_points.geometry.y, usa_points["STUSPS"]):
            if usa_points.loc[usa_points['STUSPS'] == label, 'area'].values[0] > 2.5:
                if label == 'FL':
                    x = x + 0.8
                elif label == 'LA':
                    x = x - 0.5

                if (self.nowdf.loc[self.nowdf['STUSPS'] == label, 'databin'] > (6 / 11)).all():
                    textcolor = 'w'
                else:
                    textcolor = 'k'

                texts.append(plt.text(x, y, label, fontsize=7, verticalalignment='center', fontweight=self.wgt,
                                      horizontalalignment='center', color=textcolor, fontname=self.plot_dict['myfont']))
            else:
                my_annotate(x, y, smstate[label][0], smstate[label][1], label)


    def add_basictext(self, topnum, rect2=[0.73, 0.05, 0.25, 0.35]):

        zoomscalar = 0.03
        # get the top countries

        topidx = np.flip(np.argsort(self.nowvalues))[:topnum]

        # determine if topidx contains zeros
        testconf = self.nowvalues[topidx]
        numsigdx = len(np.nonzero(testconf)[0])
        topidx = topidx[range(numsigdx)]

        ynames = list(self.nowdf.loc[topidx, 'NAME'])
        xconf = self.confdf['data'].values[topidx]
        xdeath = self.deathdf['data'].values[topidx]
        # xrec = self.recdf['data'].values[topidx]

        # hour = int(int(self.num % self.intpolnum) // (self.intpolnum / 24))
        # minute = int((self.num % (self.intpolnum / 24)) * 60 / (self.intpolnum / 24))
        # self.ax.text(0.89, 0.9, '{}:{}'.format(hour, minute), transform=self.ax.transAxes,
        #              verticalalignment='bottom', horizontalalignment='center', color='r', fontsize=12,
        #              fontname=self.plot_dict['myfont'])
        self.ax.text(0.89, 0.97, '{} {}'.format(self.day[self.num], self.month[self.num]), transform=self.ax.transAxes,
                     verticalalignment='bottom', horizontalalignment='center', color='black', fontsize=20, fontweight=self.wgt,
                     fontname=self.plot_dict['myfont'])
        self.ax.text(0.87, 0.35, 'Top {} States'.format(topnum), transform=self.ax.transAxes,
                     verticalalignment='bottom', horizontalalignment='center', color='red', fontsize=10, fontweight=self.wgt,
                     fontname=self.plot_dict['myfont'])

        #     # add flags
        nowfont = 8
        nm, fl, cs, dth = 0, 0.15, 0.5, 0.9
        ax2 = self.fig.add_axes(rect2, anchor='SE', zorder=-1)
        ax2.text((fl + nm) / 2, 0.93, 'States', transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont, fontweight=self.wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(cs, 0.93, 'Deaths(#)', transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=self.wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(dth, 0.93, 'Deaths(%)', transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=self.wgt,
                 fontname=self.plot_dict['myfont'])
        # ax2.text(rc, 0.93, 'Rcvd(%)', transform=ax2.transAxes,
        #          verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,
        #          fontname=self.plot_dict['myfont'])

        ystart = 0.76
        ystp = 0.2
        for i, cntry in enumerate(ynames):
            c = 'USstatesflag/' + list(flagname.loc[flagname['State'] == cntry, 'code'])[0] + '.png'
            offset_image(fl, ystart - i * ystp, c, ax2, zoomscalar)
            label = self.nowdf.iloc[topidx[i], 4]
            ax2.text(nm, ystart - i * ystp, '{}'.format(label), transform=ax2.transAxes,
                     verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=self.wgt,
                     fontname=self.plot_dict['myfont'])
            ax2.text(cs, ystart - i * ystp, '{}'.format(int(xdeath[i])), transform=ax2.transAxes,
                     verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=self.wgt,
                     fontname=self.plot_dict['myfont'])
            ax2.text(dth, ystart - i * ystp, '{}'.format(np.round(xdeath[i] / xconf[i] * 100, 2)),
                     transform=ax2.transAxes,
                     verticalalignment='center', horizontalalignment='center', color='red', fontsize=nowfont,fontweight=self.wgt,
                     fontname=self.plot_dict['myfont'])
            # ax2.text(rc, ystart - i * ystp, '{}'.format(np.round(xrec[i] / xconf[i] * 100, 2)), transform=ax2.transAxes,
            #          verticalalignment='center', horizontalalignment='center', color='g', fontsize=nowfont,
            #          fontname=self.plot_dict['myfont'])

        c = 'worldflag/us.png'
        offset_image(fl, ystart - topnum * 0.2, c, ax2, zoomscalar)
        ax2.text(nm, ystart - topnum * 0.2, '{}'.format('USA'), transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='black', fontsize=nowfont,fontweight=self.wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(cs, ystart - topnum * 0.2, '{}'.format(int(np.sum(self.death_all[:, self.num]))), transform=ax2.transAxes,
                 verticalalignment='center', horizontalalignment='center', color='k', fontsize=nowfont,fontweight=self.wgt,
                 fontname=self.plot_dict['myfont'])
        ax2.text(dth, ystart - topnum * 0.2,
                 '{}'.format(np.round(np.sum(self.death_all[:, self.num]) / np.sum(self.conf_all[:, self.num]), 2)),
                 transform=ax2.transAxes, verticalalignment='center', horizontalalignment='center', color='red',fontweight=self.wgt,
                 fontsize=nowfont, fontname=self.plot_dict['myfont'])
        ax2.set_axis_off()

    def add_image(self, corimg, rect1=[0, 0.7, 0.15, 0.15]):
        # inset coronavirus image
        ax1 = self.fig.add_axes(rect1, anchor='NW', zorder=-1)
        offset_image2(0.5, 0.5, corimg, ax1, 0.08, 500, 500)
        ax1.set_axis_off()

    # #     # add situational reports box
    def add_report(self, report):
        xst = 0.55
        yst = 0.-0.1
        edcolor = 'black'
        talpha = 0.7
        tcolor = 'red'
        dayidx = int(self.num // self.intpolnum) + self.plot_dict['where_to_start']
        rpnum = int(type(report.iloc[dayidx, 1]) == str)
        if rpnum == 1:
            if self.num > int(self.intpolnum * 0.1) and self.num < int(self.intpolnum * 0.9):
                text1 = report.iloc[dayidx, 1]
                if len(text1) > 60 and len(text1) <= 120:
                    tx = text1.split()
                    newtext = [" ".join(tx[: int(len(tx) / 2) + 1]), " ".join(tx[np.int(len(tx) / 2) + 1:])]
                    self.ax.text(xst, yst, '{}\n{}'.format(newtext[0], newtext[1]), transform=self.ax.transAxes,
                                 verticalalignment='center', horizontalalignment='right', color=tcolor, fontsize=7,
                                 fontname=self.plot_dict['myfont'],
                                 bbox=dict(fc='w', ec=edcolor, pad=5, alpha=talpha))
                elif len(text1) > 120:
                    tx = text1.split()
                    sptnum = np.int(len(tx) / 3)
                    newtext = [" ".join(tx[: sptnum + 1]), " ".join(tx[sptnum + 1:2 * sptnum + 1]),
                               " ".join(tx[2 * sptnum + 1:3 * sptnum + 8])]
                    self.ax.text(xst, yst, '{}\n{}\n{}'.format(newtext[0], newtext[1], newtext[2]),
                                 transform=self.ax.transAxes,
                                 verticalalignment='center', horizontalalignment='right', color=tcolor, fontsize=7,
                                 fontname=self.plot_dict['myfont'],
                                 bbox=dict(fc='w', ec=edcolor, pad=5, alpha=talpha))
                else:
                    self.ax.text(xst, yst, '{}'.format(text1), transform=self.ax.transAxes,
                                 verticalalignment='center', horizontalalignment='right', color=tcolor, fontsize=7,
                                 fontname=self.plot_dict['myfont'],
                                 bbox=dict(fc='w', ec=edcolor, pad=5, alpha=talpha))

    def saveplot(self):
        # write to video
        figname = 'G:/YoutubeChannel/CoronaVirus/CoronaImages/USA/{}{}-{}.png'.format(self.month[self.num],
                                                                                         self.day[self.num], self.numname)
        self.fig.savefig(figname, dpi=self.fig.get_dpi(), facecolor=self.fig.get_facecolor())
        plt.cla()
        self.fig.clf()
        plt.close('all')
        gc.collect()


def plotgraph(num):
    k.getmapdata(num, datatype)

    title = 'Coronavirus-USA'
    topnum = 3

    k.plotmap(title, Legend_elements, cmap, usa_points)
    k.add_basictext(topnum)
    # k.add_image(corimg)
    # k.add_report(report)
    k.saveplot()

    return f'done with {num}'



################################################################################################################
## USE PARALLEL PROCESSING TO EXPEDITE CREATION OF THOUSANDS OF FIGURES BASED ON FRAME RATE
where_to_start = 39
data_dict = dict(sort_conf=sort_conf, sort_death=sort_death, conf=conf, deat=deat)
myfont = 'Comic Sans MS'
bkgcolor = np.array([1,1,1,1])
bkgcolor[-1] = 0.95
plot_dict = dict(where_to_start=where_to_start, usa_copy=usa_copy, myfont=myfont, bkgcolor=bkgcolor, myres_w=1920,
                 myres_h=1080, mydpi=300)

k = Create_Animation(data_dict)
Create_Animation.set_variables(plot_dict)
intpolnum = 90
k.intpoldata(intpolnum)
totalnum = intpolnum * (len(sort_conf.columns[4:]) - where_to_start)
startnum = 0
num = totalnum-1
plotgraph(num)

# if __name__ == '__main__':
#     num = range(startnum, totalnum)
#
#     with concurrent.futures.ProcessPoolExecutor() as executor:
#         results = executor.map(plotgraph, num)
#
#         for result in results:
#             print(result)
