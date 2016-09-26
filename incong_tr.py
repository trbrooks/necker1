#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.80.06), Mon Aug 10 12:28:06 2015
If you publish work using this script please cite the relevant PsychoPy publications
  Peirce, JW (2007) PsychoPy - Psychophysics software in Python. Journal of Neuroscience Methods, 162(1-2), 8-13.
  Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy. Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, data, event, logging, sound, gui
from psychopy.constants import *  # things like STARTED, FINISHED
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions

# Store info about the experiment session
expName = 'incong_training'  # from the Builder filename that created this script
expInfo = {'participant':'', 'session':'001'}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False: core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName

# Setup filename for saving
filename = 'data/%s_%s_%s' %(expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=None,
    savePickle=True, saveWideText=True,
    dataFileName=filename)
#save a log file for detail verbose info
#logFile = logging.LogFile(filename+'.log', level=logging.EXP)
#logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Start Code - component code to be run before the window creation

# Setup the Window
win = visual.Window(size=(1280, 800), fullscr=True, screen=0, allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[0,0,0], rgb=[1,1,1], colorSpace='rgb',
    blendMode='avg', useFBO=True,
    )
# store frame rate of monitor if we can measure it successfully
expInfo['frameRate']=win.getActualFrameRate()
if expInfo['frameRate']!=None:
    frameDur = 1.0/round(expInfo['frameRate'])
else:
    frameDur = 1.0/60.0 # couldn't get a reliable measure so guess

# Initialize components for Routine "trial"
trialClock = core.Clock()
ISI = core.StaticPeriod(win=win, screenHz=expInfo['frameRate'], name='ISI')
hexa = visual.ImageStim(win=win, name='hexa',units='pix', 
    image='/Users/trbrooks/Desktop/Spring 2016/necker/images/hex.png', mask=None,
    ori=0, pos=[0, 0], size=[120,108],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-1.0)
cross = visual.ImageStim(win=win, name='cross',units='pix', 
    image='/Users/trbrooks/Desktop/Spring 2016/necker/images/cross.png', mask=None,
    ori=0, pos=[0, 0], size=[30,30],
    color=[1,1,1], colorSpace='rgb', opacity=0.9,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-1.0)
downy1 = visual.ImageStim(win=win, name='downy1',units=u'pix', 
    image=u'/Users/trbrooks/Desktop/Spring 2016/necker/images/downcube.png', mask=None,
    ori=0, pos=[0, 0], size=[120,108],
    color=[1,1,1], colorSpace=u'rgb', opacity=1.0,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)
downy2 = visual.ImageStim(win=win, name='downy2',units=u'pix', 
    image=u'/Users/trbrooks/Desktop/Spring 2016/necker/images/downcube.png', mask=None,
    ori=0, pos=[0, 0], size=[360,324],
    color=[1,1,1], colorSpace=u'rgb', opacity=1.0,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)
upy = visual.ImageStim(win=win, name='upy',units=u'pix', 
    image=u'/Users/trbrooks/Desktop/Spring 2016/necker/images/downcube.png', mask=None,
    ori=0, pos=[0,0], size=[120,108],
    color=[1,1,1], colorSpace=u'rgb', opacity=1.0,
    flipHoriz=True, flipVert=True,
    texRes=128, interpolate=True, depth=-3.0)
upy2 = visual.ImageStim(win=win, name='upy2',units=u'pix', 
    image=u'/Users/trbrooks/Desktop/Spring 2016/necker/images/downcube.png', mask=None,
    ori=0, pos=[0,0], size=[360,324],
    color=[1,1,1], colorSpace=u'rgb', opacity=1.0,
    flipHoriz=True, flipVert=True,
    texRes=128, interpolate=True, depth=-3.0)
mouse = event.Mouse(win=win)
x, y = [None, None]

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

#------Prepare to start Routine "trial"-------
t = 0
trialClock.reset()  # clock 
frameN = -1
routineTimer.add(150.000000)
# update component parameters for each repeat
# setup some python lists for storing info about the mouse
mouse.x = []
mouse.y = []
mouse.leftButton = []
mouse.midButton = []
mouse.rightButton = []
mouse.time = []
# keep track of which components have finished
trialComponents = []
trialComponents.append(ISI)

trialComponents.append(hexa)
trialComponents.append(downy1)
trialComponents.append(downy2)
trialComponents.append(upy2)
trialComponents.append(upy)
trialComponents.append(mouse)
for thisComponent in trialComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

#-------Start Routine "trial"-------
#create a random string of numbers from 0 to 120.
import random
temp=random.sample(range(0,160),45)
timestamp=sorted(temp)

#iterate throught the timestamps
i=iter(timestamp)
switch=i.next()

#Choose a cube to start on.
whichcube=bool(random.getrandbits(1))

#Create a list of which cube is being displayed at every time step.
cubes=[]

#At each time "temp", randomly display a stimulus, either the up cube or the down cube.


continueRoutine = True
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = trialClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    cubes.extend([whichcube])
    if trialClock.getTime() > switch:
        switch=i.next()
        whichcube= not whichcube

    # *hexa* updates
#    if t >= 0.0 and hexa.status == NOT_STARTED:
#        # keep track of start time/frame for later
#        hexa.tStart = t  # underestimates by a little under one frame
#        hexa.frameNStart = frameN  # exact frame index
#        hexa.setAutoDraw(True)
#    elif hexa.status == STARTED and t >= (0.0 + (150-win.monitorFramePeriod*0.75)): #most of one frame period left
#        hexa.setAutoDraw(False)
    
    # *downy1* updates
    if t >= 0.0 and downy1.status == NOT_STARTED:
        # keep track of start time/frame for later
        downy1.tStart = t  # underestimates by a little under one frame
        downy1.frameNStart = frameN  # exact frame index
        downy1.setAutoDraw(True)
    elif downy1.status == STARTED and t >= (0.0 + (150-win.monitorFramePeriod*0.75)): #most of one frame period left
        downy1.setAutoDraw(False)
    if downy1.status == STARTED and whichcube == True:  # only update if being drawn
        downy1.setOpacity(0.08)
    if downy1.status == STARTED and whichcube == False:  # only update if being drawn
        downy1.setOpacity(1)

    # *upy* updates
    if t >= 0.0 and upy.status == NOT_STARTED:
        # keep track of start time/frame for later
        upy.tStart = t  # underestimates by a little under one frame
        upy.frameNStart = frameN  # exact frame index
        upy.setAutoDraw(True)
    elif upy.status == STARTED and t >= (0.0 + (150-win.monitorFramePeriod*0.75)): #most of one frame period left
        upy.setAutoDraw(False)
    if upy.status == STARTED and whichcube == False:  # only update if being drawn
        upy.setOpacity(0.08)
    if upy.status == STARTED and whichcube == True:  # only update if being drawn
        upy.setOpacity(1)

  # *cross* updates
    if t >= 0.0 and cross.status == NOT_STARTED:
        # keep track of start time/frame for later
        cross.tStart = t  # underestimates by a little under one frame
        cross.frameNStart = frameN  # exact frame index
        cross.setAutoDraw(True)
    elif cross.status == STARTED and t >= (0.0 + (150-win.monitorFramePeriod*0.75)): #most of one frame period left
        cross.setAutoDraw(False)


# *mouse* updates
    if t >= 0.0 and mouse.status == NOT_STARTED:
        # keep track of start time/frame for later
        mouse.tStart = t  # underestimates by a little under one frame
        mouse.frameNStart = frameN  # exact frame index
        mouse.status = STARTED
        event.mouseButtons = [0, 0, 0]  # reset mouse buttons to be 'up'
    elif mouse.status == STARTED and t >= (0.0 + (150.0-win.monitorFramePeriod*0.75)): #most of one frame period left
        mouse.status = STOPPED
    if mouse.status == STARTED:  # only update if started and not stopped!
        buttons = mouse.getPressed()
        x, y = mouse.getPos()
        mouse.x.append(x)
        mouse.y.append(y)
        mouse.leftButton.append(buttons[0])
        mouse.midButton.append(buttons[1])
        mouse.rightButton.append(buttons[2])
        mouse.time.append(trialClock.getTime())
        
    # *downy2* updates
    if t >= 0.0 and downy2.status == NOT_STARTED:
        # keep track of start time/frame for later
        downy2.tStart = t  # underestimates by a little under one frame
        downy2.frameNStart = frameN  # exact frame index
        downy2.setAutoDraw(True)
    elif downy2.status == STARTED and t >= (0.0 + (150-win.monitorFramePeriod*0.75)): #most of one frame period left
        downy2.setAutoDraw(False)
    if downy2.status == STARTED and mouse.leftButton[-1] != 0:  # only update if being drawn
        downy2.setOpacity(0.08)
    if downy2.status == STARTED and mouse.leftButton[-1] == 0:  # only update if being drawn
        downy2.setOpacity(0.5)
     
    # *upy2* updates
    if t >= 0.0 and upy2.status == NOT_STARTED:
    # keep track of start time/frame for later
        upy2.tStart = t  # underestimates by a little under one frame
        upy2.frameNStart = frameN  # exact frame index
        upy2.setAutoDraw(True)
    elif downy2.status == STARTED and t >= (0.0 + (150-win.monitorFramePeriod*0.75)): #most of one frame period left
        upy2.setAutoDraw(False)
    if upy2.status == STARTED and mouse.leftButton[-1] == 0:  # only update if being drawn
        upy2.setOpacity(0.08)
    if upy2.status == STARTED and mouse.leftButton[-1] != 0:  # only update if being drawn
        upy2.setOpacity(0.5)


# *ISI* period
    if t >= 0.0 and ISI.status == NOT_STARTED:
        # keep track of start time/frame for later
        ISI.tStart = t  # underestimates by a little under one frame
        ISI.frameNStart = frameN  # exact frame index
        ISI.start(0.5)
    elif ISI.status == STARTED: #one frame should pass before updating params and completing
        ISI.complete() #finish the static period
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        routineTimer.reset()  # if we abort early the non-slip timer needs reset
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in trialComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

#-------Ending Routine "trial"-------
for thisComponent in trialComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# store data for thisExp (ExperimentHandler)
thisExp.addData('mouse.leftButton', mouse.leftButton)
thisExp.addData('mouse.time', mouse.time)
thisExp.addData('time.stamps',timestamp)
thisExp.addData('cubes',cubes)
thisExp.addData('iscong', 0)
thisExp.addData('istest', 0)
thisExp.nextEntry()
win.close()
core.quit()
