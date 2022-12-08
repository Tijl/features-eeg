#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from psychopy import core, event, visual, gui
import os,random,sys,math,json,requests
from glob import glob
import pandas as pd
import numpy as np

# debug things
debug_testsubject = 1
debug_usedummytriggers = 1
debug_windowedmode = 1
debug_save_screenshots = 0

if debug_testsubject:
    subjectnr = 0
else:
    # Get subject info
    subject_info = {'Subject number':''}
    if not gui.DlgFromDict(subject_info,title='Enter subject info:').OK:
        print('User hit cancel at subject information')
        exit()
    try:
        subjectnr = int(subject_info['Subject number'])
    except:
        raise

outfn = 'sub-%02i_task-rsvp_events.csv'%subjectnr
if not debug_testsubject and os.path.exists(outfn):
    raise Exception('%s exists'%outfn)

if not debug_usedummytriggers:
    from psychopy import parallel

random.seed(subjectnr)

nsequence = 40 #repeats per stimulus

refreshrate = 60
fixationduration = 1 #how long to pause before starting the sequence (s)
stimdurations = [.01, .01]  #stimulus duration 
soadurations = [.15 , .05]  #interstimulus time

trigger_stimon = 1
trigger_stimoff = 2
trigger_sequencestart = 3
trigger_duration = 0.010  
trigger_port = 0xdff8

with open('webhook_url') as w:
    webhook_url=w.readline()

#feature levels
feature_ori = [22.5,67.5,112.5,157.5]
feature_sf = [.055, .04, .025, .01]
feature_color = [
	( 66, 10,104),
	(147, 38,103),
	(221, 81, 58),
	(252,165, 10),
    ]
feature_contrast = [0.9, 0.7, 0.5, 0.3]

stims = np.array(np.meshgrid(range(4), range(4), range(4), range(4))).T.reshape(-1,4)

nstimuli = len(stims)
stimnumber = list(range(nstimuli))

allsequence = []
for s in range(nsequence):
    allsequence += random.sample(stimnumber,nstimuli)
allsequence += allsequence #double for the two speed conditions

eventlist = pd.DataFrame(allsequence,columns=['stimnumber'])
eventlist['f_ori'] = [stims[x,0] for x in eventlist['stimnumber']]
eventlist['f_sf'] = [stims[x,1] for x in eventlist['stimnumber']]
eventlist['f_color'] = [stims[x,2] for x in eventlist['stimnumber']]
eventlist['f_contrast'] = [stims[x,3] for x in eventlist['stimnumber']]
eventlist['feature_ori'] = [feature_ori[x] for x in eventlist['f_ori']]
eventlist['feature_sf'] = [feature_sf[x] for x in eventlist['f_sf']]
eventlist['feature_color'] = [feature_color[x] for x in eventlist['f_color']]
eventlist['feature_contrast'] = [feature_contrast[x] for x in eventlist['f_contrast']]
eventlist['blocksequencenumber'] = [math.floor(x/nstimuli%nsequence) for x in range(len(eventlist))]
eventlist['sequencenumber'] = [math.floor(x/(nstimuli)) for x in range(len(eventlist))]
eventlist['speedcondition'] = [math.floor(x/nstimuli/nsequence) for x in range(len(eventlist))]
eventlist['stimduration'] = [stimdurations[x] for x in eventlist['speedcondition']]
eventlist['soaduration'] = [soadurations[x] for x in eventlist['speedcondition']]

nsequencetotal = eventlist['sequencenumber'].iloc[-1]+1

#add targets
istarget = []
for s in range(nsequencetotal):
    nr = len(eventlist[eventlist['sequencenumber']==s])
    ntargets = random.randint(2,4) #2 to 4 targets in a seq. 
    tt=[False for x in range(nr)]
    targetpos=[1, 1]
    while any(np.diff(targetpos)<20):
        targetpos = sorted(random.sample(range(10,nr-10),ntargets))
    for p in targetpos:
        tt[p]=True
    istarget+=tt
eventlist['istarget']=istarget

#concat, then resort by index to insert target stims
# etarget = eventlist.loc[istarget,:]
# etarget.loc[:,'istarget']=True
# eventlist = pd.concat((eventlist,etarget)).sort_index(kind='mergesort').reset_index(drop=True)

# shuffle sequences
sorder = random.sample(range(nsequencetotal),nsequencetotal)
eventlist['sequencenumber'] = [sorder[x] for x in eventlist['sequencenumber']]
eventlist = eventlist.sort_values('sequencenumber',kind='mergesort').reset_index(drop=True)
# add presentationnumber in sequence
eventlist['presentationnumber'] = [1 for x in range(len(eventlist))]
eventlist['presentationnumber'] = eventlist.groupby(['sequencenumber'])['presentationnumber'].cumsum()-1

def writeout(eventlist):
    with open(outfn,'w') as out:
        eventlist.to_csv(out,index_label='eventnumber')

writeout(eventlist)

# =============================================================================
# %% START
# =============================================================================
try:
    if debug_windowedmode:
        win=visual.Window([800,800],units='pix')
    else:
        win=visual.Window(units='pix',fullscr=True)
    mouse = event.Mouse(visible=False)

    fixation = visual.ImageStim(win, 'eeg_fixation.png', size=24,
         name='fixation', autoLog=False)
    fixationtarget = visual.ImageStim(win, 'eeg_fixation_target.png', size=24,
         name='fixationtarget', autoLog=False)
    
    querytext = visual.TextStim(win,text='',pos=(0,200),name='querytext')
    progresstext = visual.TextStim(win,text='',pos=(0,100),name='progresstext')
    sequencestarttext = visual.TextStim(win,text='',pos=(0,50),name='sequencestarttext')
    
    filesep='/'
    if sys.platform == 'win32':
        filesep='\\'
        
    def check_abort(k):
        if k and k[0][0]=='q':
            raise Exception('User pressed q')
            
    screenshotnr = 0
    def take_screenshot(win):
        global screenshotnr 
        screenshotnr += 1
        win.getMovieFrame()
        win.saveMovieFrames('screenshots/screen_%05i.png'%screenshotnr)
        
    def make_stimtex():
        stimtex=[]
        for i in range(nstimuli):
            s = visual.GratingStim(win, tex='sin', mask='circle', size=256, colorSpace='rgb255',
                ori=feature_ori[stims[i,0]],
                sf=feature_sf[stims[i,1]],
                color=feature_color[stims[i,2]],
                contrast=feature_contrast[stims[i,3]],
                phase=random.random()
                )
            s.name = 'stim%03i_ori%i_sf%i_color%i_contrast%i.png'%(i,stims[i,0],stims[i,1],stims[i,2],stims[i,3])
            stimtex.append(s)
        return stimtex

    if debug_save_screenshots:
        stimtex = make_stimtex()
        for s in stimtex:
            s.draw()
            win.flip()
            win.getMovieFrame()
            print(s.name)
            print(s)
            win.saveMovieFrames('stimuli'+filesep+s.name)
        
    def send_dummy_trigger(trigger_value):
        core.wait(trigger_duration)
            
    def send_real_trigger(trigger_value):
        trigger_port.setData(trigger_value)
        core.wait(trigger_duration)
        trigger_port.setData(0)
    
    if debug_usedummytriggers:
        sendtrigger = send_dummy_trigger
    else:
        trigger_port = parallel.ParallelPort(address=trigger_port)
        trigger_port.setData(0)
        sendtrigger = send_real_trigger

    nevents = len(eventlist)
    sequencenumber = -1
    for eventnr in range(nevents):
        first = eventlist['sequencenumber'].iloc[eventnr]>sequencenumber
        if first: #start of sequence
            writeout(eventlist)
            
            sequencenumber = eventlist['sequencenumber'].iloc[eventnr]
            
            stimduration = eventlist['stimduration'].iloc[eventnr] - .5/refreshrate
            soaduration = eventlist['soaduration'].iloc[eventnr] - .5/refreshrate
            
            last_target = -99
            correct=0
            
            fixation.draw()
            win.flip()

            stimtex = make_stimtex()
                
            if sequencenumber:
                idx = eventlist['sequencenumber']==(sequencenumber-1)
                nt = sum(eventlist['istarget'][idx])
                nc = sum(eventlist['correct'][idx])
                slack_data={'text':'sub-%02i seq %i/%i (V1) hits %i/%i <@tijlgrootswagers><@U9C24ECQ7>'%(subjectnr,sequencenumber,nsequencetotal,nc,nt),'channel':'#eeglab','username':'python'}
                if debug_testsubject:
                    print(slack_data)
                else:
                    try:
                        response = requests.post(webhook_url, data=json.dumps(slack_data),headers={'Content-Type': 'application/json'})
                    except:
                        pass
                
            progresstext.text = '%i / %i'%(1+sequencenumber,nsequencetotal)
            progresstext.draw()
            sequencestarttext.text = 'Press leftmost key to start the sequence'
            sequencestarttext.draw()
            fixation.draw()
            win.flip()
            k=event.waitKeys(keyList='afq', modifiers=False, timeStamped=True) #the keys of afq correspond to outer two buttons on the keypad/ start or quit
            check_abort(k)
            fixation.draw()
            time_fixon = win.flip()
            sendtrigger(trigger_sequencestart)
            while core.getTime() < time_fixon + fixationduration - .5/refreshrate:pass
        
        response=0
        rt=0
        istarget = eventlist['istarget'].iloc[eventnr]
        stimnum = eventlist['stimnumber'].iloc[eventnr]
        stim = stimtex[stimnum]
        stim.draw()
        if istarget:
            fixationtarget.draw()
        else:
            fixation.draw()
        time_stimon=win.flip()
        sendtrigger(trigger_stimon)
        
        if istarget:
            last_target=time_stimon
        
        while core.getTime() < time_stimon + stimduration:pass
        if istarget:
            fixationtarget.draw()
        else:
            fixation.draw()
        time_stimoff=win.flip()
        sendtrigger(trigger_stimoff)
        
        #get response
        correct=0
        k=event.getKeys(keyList='sdq', modifiers=False, timeStamped=True)   #keys of sdq for target selection correspond to two inner keys on keypad/ target selection
        if k:
            check_abort(k)
            response=1
            rt=k[0][1]
            correct = rt-last_target < 1.5
        
        eventlist.at[eventnr, 'feature_phase'] = stim.phase[0]
        eventlist.at[eventnr, 'stimname'] = stim.name
        eventlist.at[eventnr, 'response'] = int(response)
        eventlist.at[eventnr, 'rt'] = rt-last_target if correct else 0
        eventlist.at[eventnr, 'correct'] = int(correct)
        eventlist.at[eventnr, 'time_stimon'] = time_stimon
        eventlist.at[eventnr, 'time_stimoff'] = time_stimoff
        eventlist.at[eventnr, 'stimdur'] = time_stimoff-time_stimon
        while core.getTime() < time_stimon + soaduration:pass

finally:
    writeout(eventlist)
    sequencestarttext.text='Experiment finished!'
    sequencestarttext.draw()
    win.flip()
    core.wait(1)
    win.close()
    exit()

