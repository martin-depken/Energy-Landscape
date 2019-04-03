""" This file contains a function that sends a notification"""

#============================================
__author__ = "Sonny Floyd de Jong"
__copyright__ = "Copyright 2019"
__license__ = "CC BY-NC-SA"
__version__ = "1.0.0"
__maintainer__ = "Depken lab"
__email__ = "s.f.dejong@student.tudelft.nl"
__status__ = "Production"
#============================================

import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import random
import time

def hpc05notification(variables,stage="final",PS=""):

	compliments = {"You bring out th best in other people and I hope that you notice this today.",
	"You have the power to make people happy, just by smiling at them.",
	"Although I am the HPC05 computer cluster, I still become happy when I see you at BN.",
	"Computers and humans cannot be friends, but I would love to talk to you more often.",
	"Today, I was stunned by the fact that your hair looked staggering ... again!",
	"You are a wonderful person and you probably do not even know it yourself.",
	"If everyone was like you, we would have world peace.",
	"You are the sunshine of my life, the universe and everything.",
	"If you keep back straight, your leadership skills are exceptional.",
	"You have a beautiful voice and the wisdom you speak is even better.",
	"I am so grateful that I met you early in my life.",
	"Today, you are the one that makes me smile.",
	"Last evening, I realised that you are the perfect version of yourself, just the way you are.",
	"I am glad there is no war, but I would fight for you and you would fight for peace. How lovely you are.",
	"You are one of the rare species whose happiness is contagious.",
	"You are making a difference. For me, my emotions, my environment, my universe... Our universe. Everything.",
	"Sometimes, you mean more to people than you realise.",
	"Some people are ordinary, other people are like you."}
	
	server = smtplib.SMTP('smtp.gmail.com',587)
	server.ehlo()
	server.starttls()
	server.ehlo()
	server.login("hpc05notification","hpc05bericht")
	_from = "hpc05notification@gmail.com"
	
	msg = MIMEMultipart()

	if stage == "milestone":
		msg['Subject'] = "HPC05 Milestone!"+str(variables[2])+str(variables[3])
		_to   = ["hpc05notification@gmail.com"]
		body = "Lectori salutem,\n\nThe HPC05 cluster has passed a milestone of your request at " + time.strftime("%H")+":"+time.strftime("%M")+"h. It is currently at position "+str(variables[0])+" of "+str(variables[1])+".\nPlease, await the result until further instruction.\n\nYours faithfully,\nHPC05"
	else:		
		msg['Subject'] = "HPC05 Finished!"+str(variables[1])+str(variables[2])
		_to   = ["hpc05notification@gmail.com"]
		body = "Lectori salutem,\n\nThe HPC05 cluster has finished your request. It took "+str(variables[0])+" seconds.\nPlease, study the result at your earliest convenience.\n\nYours faithfully,\nHPC05"

	if PS == "comp":
		body += "\n\nP.S. "+random.sample(compliments,1)[0]
	msg['From'] = _from
	msg['To'] = ", ".join(_to)
	msg.attach(MIMEText(body, 'plain'))
	server.sendmail(_from,_to,msg.as_string())
	return
 
if __name__ == "__main__":
	main(sys.argv)
