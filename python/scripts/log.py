import os
import logging
import logging.config

HOST = 'smtp.cruk.cam.ac.uk'
FROM = 'anne.pajon@cruk.cam.ac.uk'
TO = 'anne.pajon@cruk.cam.ac.uk'
SUBJECT = '[ERROR] Genome Editing Project'

# logging definition
LOGGING = {
    'version': 1,
    'disable_existing_loggers': True,
    'formatters': {
        'verbose': {
            'format': '%(asctime)s %(name)-24s %(levelname)-8s: %(message)s'
        },
        'simple': {
            'format': '%(name)-24s %(levelname)-8s: %(message)s'
        },
    },
    'handlers': {
        'console': {
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple'
        },
        'info_file': {
            'level': 'DEBUG',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': 'info.log',
            'maxBytes': 1000000,
            'backupCount': 5,
            'formatter': 'verbose',
        },
        'error_file': {
            'level': 'ERROR',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': 'errors.log',
            'maxBytes': 1000000,
            'backupCount': 5,
            'formatter': 'verbose',
        },
        'email': {
            'level': 'ERROR',
            'class': 'logging.handlers.SMTPHandler',
            'mailhost': HOST,
            'fromaddr': FROM,
            'toaddrs': [TO],
            'subject': SUBJECT,
            'formatter': 'verbose',
        },
    },
    'loggers': {
        'dnascissors': {
            'handlers': ['console', 'info_file', 'error_file'],
            'propagate': True,
            'level': 'DEBUG',
        },
    }
}


def get_custom_logger(logfile=None, noemail=False):
    if logfile:
        if not os.path.exists(os.path.dirname(logfile)):
            os.makedirs(os.path.dirname(logfile))
        LOGGING['handlers']['info_file']['filename'] = logfile
        LOGGING['handlers']['error_file']['filename'] = logfile + ".errors"
    if not noemail:
        LOGGING['loggers']['dnascissors']['handlers'].append('email')
    logging.config.dictConfig(LOGGING)
    return logging.getLogger('dnascissors')
