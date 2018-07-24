import smtplib
import sys

def emailalert(message,password,email):
    # SMTP_SSL Example
    server_ssl = smtplib.SMTP_SSL("smtp.gmail.com", 465)
    server_ssl.ehlo() # optional, called by login()
    server_ssl.login(email, password)  
    # ssl server doesn't support or need tls, so don't call server_ssl.starttls() 
    server_ssl.sendmail(email, email, message)
    #server_ssl.quit()
    server_ssl.close()
    print('successfully sent the mail')


if __name__ == "__main__":
    message=sys.argv[1]
    password=sys.argv[2]
    email=sys.argv[3]
    emailalert(message,password,email)