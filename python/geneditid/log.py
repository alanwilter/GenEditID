import os
import logging
import logging.config

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
    },
    'loggers': {
        'geneditid': {
            'handlers': ['console', 'info_file'],
            'propagate': True,
            'level': 'DEBUG',
        },
    }
}


def get_custom_logger(logfile=None):
    if logfile:
        logdir = os.path.dirname(logfile)
        if logdir:
            if not os.path.exists(logdir):
                os.makedirs(logdir)
        LOGGING['handlers']['info_file']['filename'] = logfile
    logging.config.dictConfig(LOGGING)
    return logging.getLogger('geneditid')
