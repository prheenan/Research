#!/usr/bin/env python

'''
COMMAND MACRO PLUGIN FOR HOOKE

Records, saves and executes batches of commands
(c)Alberto Gomez-Casado 2008
'''

import libhookecurve as lhc
import libinput as linput
import os.path
import string

class macroCommands:

	currentmacro=[]
	pause=0
	auxprompt=[]
	macrodir=[]
	

	def _plug_init(self):
		self.currentmacro=[]
		self.auxprompt=self.prompt
		self.macrodir=self.config['workdir']
		if not os.path.exists(os.path.join(self.macrodir,'macros')):
                    try:
                        os.mkdir('macros')
                    except:
                        print 'Warning: cannot create macros folder.'
                        print 'Probably you do not have permissions in your Hooke folder, use macro at your own risk.'
                self.macrodir=os.path.join(self.macrodir,'macros')

	def collect(self):
				
		print 'Enter STOP / PAUSE to go back to normal mode\nUNDO to remove last command'
		line=[]
		while not(line=='STOP' or line=='PAUSE'):
			line=raw_input('hooke (macroREC): ')
			if line=='PAUSE':
				self.pause=1
				self.prompt='hooke (macroPAUSE): '
				break
			if line=='STOP':
				self.prompt=self.auxprompt
				self.do_recordmacro('stop')
				break
			if line=='UNDO':
				self.currentmacro.pop()
				continue
			param=line.split()

			#FIXME check if accessing param[2] when it doesnt exist breaks something
			if param[0] =='export':
				exportline=param[0]+' __curve__ '
				if len(param)==3:
					exportline=exportline+param[2]
				self.currentmacro.append(exportline)
				self.onecmd(line)
				continue
			
			if param[0] =='txt':
				exportline=param[0]
				if len(param)==3:
					exportline=exportline+' '+param[2]
				exportline=exportline+'__curve__'
				self.currentmacro.append(exportline)
				self.onecmd(line)
				continue

			self.onecmd(line)
			
			self.currentmacro.append(line)
		

	def do_recordmacro(self, args):
		'''RECORDMACRO
		Stores input commands to create script files
		-------
		Syntax: recordmacro [start / stop]
		If a macro is currently paused start resumes recording
		'''
		
		
		if len(args)==0:
			args='start'

		if args=='stop':
			self.pause=0
			self.prompt=self.auxprompt
			if len(self.currentmacro) != 0: 
				answer=linput.safeinput('Do you want to save this macro? ',['y'])
				if answer[0].lower() == 'y':
					self.do_savemacro('')
				else:
					print 'Macro discarded'
					self.currentmacro=[]
			else:
				print 'Macro was empty'	

		if args=='start':	

			if self.pause==1:
				self.pause=0	
				self.collect()	
			else:
				if len(self.currentmacro) != 0: 
					answer=linput.safeinput('Another macro is already beign recorded\nDo you want to save it?',['y'])
					if answer[0].lower() == 'y':
						self.do_savemacro('')
					else:
						print 'Old macro discarded, you can start recording the new one'
			
				self.currentmacro=[]
				self.collect()
		

	def do_savemacro(self, macroname):

		'''SAVEMACRO
		Saves previously recorded macro into a script file for future use
		-------
		Syntax: savemacro [macroname]
		If no macroname is supplied one will be interactively asked 
		'''

		saved_ok=0
		if self.currentmacro==None:
			print 'No macro is being recorded!'
			return 0
		if len(macroname)==0: 
			macroname=linput.safeinput('Enter new macro name: ')
			if len(macroname) == 0:
				print 'Invalid name'
				
		macroname=os.path.join(self.macrodir,macroname+'.hkm')
		if os.path.exists(macroname):
			overwrite=linput.safeinput('That name is in use, overwrite?',['n'])
			if overwrite[0].lower()!='y':
				print 'Cancelled save'
				return 0
		txtfile=open(macroname,'w+')
		self.currentmacro='\n'.join(self.currentmacro)
		txtfile.write(self.currentmacro)
		txtfile.close()
		print 'Saved on '+macroname
		self.currentmacro=[]

	def do_execmacro (self, args):
		
		'''EXECMACRO
		Loads a macro and executes it over current curve / playlist
		-----
		Syntax: execmacro macroname [playlist] [v]

		macroname.hkm should be present at [hooke]/macros directory
		By default the macro will be executed over current curve
		passing 'playlist' word as second argument executes macroname 
		over all curves
		By default curve(s) will be processed silently, passing 'v' 
		as second/third argument will print each command that is 
		executed

		Note that macros applied to playlists should end by export
		commands so the processed curves are not lost
		'''
		verbose=0
		cycle=0
		curve=None		

		if len(self.currentmacro) != 0:
			print 'Warning!: you are calling a macro while recording other'
		if len(args) == 0:
			print 'You must provide a macro name'
			return 0
		args=args.split()

		#print 'args ' + ' '.join(args)
		
		if len(args)>1:
			if args[1] == 'playlist':
				cycle=1
				print 'Remember! macros applied over playlists should include export orders'
				if len(args)>2 and args[2] == 'v':
					verbose=1
			else:
				if args[1] == 'v':
					verbose=1	
		#print cycle
		#print verbose	

		macropath=os.path.join(self.macrodir,args[0]+'.hkm')
		if not os.path.exists(macropath):
			print 'Could not find a macro named '+macropath
			return 0
		txtfile=open(macropath)
		if cycle ==1: 
			#print self.current_list
			for item in self.current_list:
				self.current=item
				self.do_plot(0)        

				for command in txtfile:

					if verbose==1:
						print 'Executing command '+command
					testcmd=command.split()
					w=0
					for word in testcmd:
						if word=='__curve__':
							testcmd[w]=os.path.splitext(os.path.basename(item.path))[0]
						w=w+1
					self.onecmd(' '.join(testcmd))
				self.current.curve.close_all()
				txtfile.seek(0)
		else:
			for command in txtfile:
					testcmd=command.split()
					w=0
					for word in testcmd:
						if word=='__curve__':
							w=w+1
							testcmd[w]=os.path.splitext(os.path.basename(self.current.path))[0]+'-'+string.lstrip(os.path.splitext(os.path.basename(self.current.path))[1],'.')
					if verbose==1:
						print 'Executing command '+' '.join(testcmd)
					self.onecmd(' '.join(testcmd))
		
	




