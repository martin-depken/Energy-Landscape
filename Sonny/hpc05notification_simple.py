""" This file contains a function that sends a notification from the HPC05 cluster."""

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

def hpc05notification():
	"""
	1. Add the desired email addresses in the variable _to. It can be multiple or one.
	2. Import the only function of this file into Python, for instance by the command:
	   from hpc05notification_simple import hpc05notification
	3. Execute hpc05notification() in your code.
	
	"""

	# FILL IN YOUR RECIPIENTS:
	_to   = ["s.f.dejong@student.tudelft.nl", "hpc05notification@gmail.com"]
	
	
	server = smtplib.SMTP('smtp.gmail.com',587)
	server.ehlo()
	server.starttls()
	server.ehlo()
	server.login("hpc05notification","hpc05bericht")
	_from = "hpc05notification@gmail.com"
	msg = MIMEMultipart()
	msg['Subject'] = "Your HPC05 Finished!"
	body = "Lectori salutem,\n\nThe HPC05 cluster has finished your request.\nPlease, study the result at your earliest convenience.\n\nYours faithfully,\nHPC05"
	msg['From'] = _from
	msg['To'] = ", ".join(_to)
	msg.attach(MIMEText(body, 'plain'))
	server.sendmail(_from,_to,msg.as_string())
	return
 
if __name__ == "__main__":
	main(sys.argv)
