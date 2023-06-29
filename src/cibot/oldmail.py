"""
Old routine for sending messages through smtplib, using password authentication.
This no longer works due to gmail's 2022 API changes. Thanks google.

@author: Nick Gibbons
"""
import smtplib

TEMPLATE = """\
From: {sender}
To: {to}
Subject: {subject}

{body}

I am a bot and this action was performed automatically!
Please contact {handler} if you have any questions or concerns.
"""

class Emailer(object):
    def __init__(self, username, passwordfile):
        self.username = username
        with open(passwordfile) as fp:
            self.password = fp.read()

    def send(self, to, message):
        try:
            server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
            server.ehlo()
            server.login(self.username, self.password)
            server.sendmail(self.username, to, message)
            server.close()
            print("Sent!")
        except:
            print("Failed to send mail.")


if __name__=='__main__':
    username = 'eilmer4bot@gmail.com'
    passwordfile = '/home/eilmer/.gmail_password'

    to = ['n.gibbons@uq.edu.au']
    subject = 'Commit Test!'
    body = """\
commit 5f0340f25249728c8d0239e2a9cf8a8d31679fa0 (HEAD -> master)
Author: Nick N. Gibbons <n.gibbons@uq.edu.au>
Date:   Mon Nov 9 17:59:21 2020 +1000

    Added some basic mail sending
    """

    message = TEMPLATE.format(sender=username, to=','.join(to), subject=subject, body=body, handler=to[0])
    print(message)
    emailer = Email(username, passwordfile)
    emailer.send(to, message)


