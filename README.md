# ANIMATED TIME SERIES VISUALIZATION OF COVID-19 PATTERNS ACROSS THE WORLD

## Overview
The COVID-19 pandemic, caused by the novel coronavirus, has affected people in almost every country in the world, making it a truly global phenomenon. The virus is highly contagious, is able to spread rapidly and can cause serious illness and death, making it a major public health concern.The pandemic has had major impacts on the global economy, leading to widespread job losses and business closures. It has also forced people to make major changes to their daily routines, including staying at home and avoiding large gatherings, which has been a significant shift in people's lives.

In this project, I used solely python to conduct data cleaning and exploration and created an animated time series visualization of the spread of the virus across countries during the pandemic phase in 2020. In python, I used several techniques and libraries including parallel processing to create tens of thousands of HD images, matloplib and geopandas for beautiful edited images.  I shared most of my work on my youtube channel here https://www.youtube.com/@statsmania504 getting 86 subscribers and thousands of views. 

This repo contains python script to create data visualization for COVID data

There are two types of visualizations created here;

__1. Bar graphs:__ This type of visualization uses bar graphs to track how the top N countries/states based on a particular epidemiological metric (e.g. death rate) of COVID-19 change during a period of time.

Examples of the visualizations were shared in youtube videos here:

a. Top 20 Countries by COVID-19 Deaths per 1 million population (Jan 22nd to Jun 9th, 2020)
https://www.youtube.com/watch?v=y4irhI3GCXo&t=20s&ab_channel=StatsMania

b. Top 20 US States by Total COVID-19 Cases per person and per 1 million population (Feb 15 - Jun 11)
https://www.youtube.com/watch?v=V34ELCXAPGA&ab_channel=StatsMania



__2. Maps:__ This type of visualization uses geomaps (from geopandas python library) to track how a particular epidemiological metric (e.g. death rate) of COVID-19 affected a particular geographical area during a period of time.

Examples of the visualizations were shared in youtube videos here:

a. Map timelapse of total COVID-19 Cases in the US from March 1st to June 12th
https://www.youtube.com/watch?v=4My3Fs_ulTw&t=33s&ab_channel=StatsMania

b. Global Map Timelapse of Active COVID-19 Cases from Feb 1st to June 9th, 2020.
https://www.youtube.com/watch?v=Wg6ih6t2JNc&ab_channel=StatsMania
  
## Code

There are two steps here:
1. The first step involves create several thousands of figures based on the duration and frame rate of the video
2. The second involves involves combining the generated figures into a video.
N.B. While most of the details in the video were created using python, a few aspects such as the introduction frames were created using video editing softwares

___For Bar Graps:___

Step 1: CoronavirusBarActive-World.py, CoronaVirusBar-PerCapita-World.py, CoronaVirusBar-Raw+PerCapita_World.py

Step 2: makingvideo.py

___For Maps:___

Step 1:CoronavirusMapActive-World.py, CoronavirusMapActive-USA.py

Step 2:makingmap.py


