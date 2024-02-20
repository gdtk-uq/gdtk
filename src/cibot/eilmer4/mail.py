"""
Gmail sending script with much help from stack overflow:
 - https://stackoverflow.com/questions/37201250/sending-email-via-gmail-python

Requires:
 pip3 install google-api-python-client oauth2client

Notes:
 Invoke this file using $ python3 mail.py to bring up a sign in window to
 Google. At that point you can log into the eilmer4bot account and click agree
 to allow mail sending. This will fresh the bot's token for one week.

@author: Nick Gibbons
"""
import httplib2
import os
from oauth2client import client, tools, file
import base64
from email import encoders
from email.mime.text import MIMEText
from apiclient import errors, discovery

TEMPLATE = """\

{body}

I am a bot and this action was performed automatically!
Please contact {handler} if you have any questions or concerns.
"""

def get_credentials():
    # If needed create folder for credential
    home_dir = os.path.expanduser('~')
    credential_dir = os.path.join(home_dir, '.credentials')
    if not os.path.exists(credential_dir):
        os.makedirs(credential_dir)
    credential_path = os.path.join(credential_dir, 'cred_send_mail.json')

    store = file.Storage(credential_path)
    credentials = store.get()

    if not credentials or credentials.invalid:
        CLIENT_SECRET_FILE = '/home/testbot/client_secret.json'
        APPLICATION_NAME = 'Eilmer Continuous Integration'
        SCOPES = 'https://www.googleapis.com/auth/gmail.send'

        # Create a flow object. (it assists with OAuth 2.0 steps to get user authorization + credentials)
        flow = client.flow_from_clientsecrets(CLIENT_SECRET_FILE, SCOPES)
        flow.user_agent = APPLICATION_NAME
        credentials = tools.run_flow(flow, store)
    return credentials

## Get creds, prepare message and send it
def create_message_and_send(sender, to, subject,  message_text_plain):
    print("Getting credentials...")
    credentials = get_credentials()

    print("Authorising credentials...")
    # Create an httplib2.Http object to handle our HTTP requests, and authorize it using credentials.authorize()
    http = httplib2.Http()
    # http is the authorized httplib2.Http() 
    http = credentials.authorize(http)        #or: http = credentials.authorize(httplib2.Http())
    print("Getting service...")
    service = discovery.build('gmail', 'v1', http=http)

    message_without_attachment = create_message_without_attachment(sender, to, subject, message_text_plain)
    print("Sending Message...")
    send_Message_without_attachment(service, "me", message_without_attachment, message_text_plain)
    print("Sent")

def create_message_without_attachment (sender, to, subject, message_text_plain):
    #Create message container
    message = MIMEText(message_text_plain, 'plain')
    message['Subject'] = subject
    message['From'] = sender
    message['To'] = to

    raw_message_no_attachment = base64.urlsafe_b64encode(message.as_bytes())
    raw_message_no_attachment = raw_message_no_attachment.decode()
    body  = {'raw': raw_message_no_attachment}
    return body

def send_Message_without_attachment(service, user_id, body, message_text_plain):
    try:
        message_sent = (service.users().messages().send(userId=user_id, body=body).execute())
        message_id = message_sent['id']
        print (f'Message sent (without attachment) \n\n Message Id: {message_id}\n\n Message:\n\n {message_text_plain}')
    except errors.HttpError as error:
        print (f'An error occurred: {error}')

if __name__ == '__main__':
    to = ["n.gibbons@uq.edu.au", "kyle.damm@uqconnect.edu.au"]
    sender = "eilmer4bot@gmail.com"
    subject = "Commit Test!"
    body = """\
commit 5f0340f25249728c8d0239e2a9cf8a8d31679fa0 (HEAD -> master)
Author: Nick N. Gibbons <n.gibbons@uq.edu.au>
Date:   Mon Nov 9 17:59:21 2020 +1000

    Added some basic mail sending
"""
    to = ','.join(to)
    create_message_and_send(sender, to, subject, body)
