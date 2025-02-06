# logger.py

import logging
import sys

# Configure the base logger
logger = logging.getLogger("orthoxml")
logger.setLevel(logging.DEBUG)

# Create formatters
console_formatter = logging.Formatter(
    fmt='%(asctime)s - %(module)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# Console handler
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(console_formatter)

# Add handlers
logger.addHandler(console_handler)

def get_logger(name: str = None) -> logging.Logger:
    """Get a logger instance.
    
    Args:
        name: Optional name for the logger. If None, returns root logger.
        
    Returns:
        Logger instance
    """
    if name:
        return logging.getLogger(f"orthoxml.{name}")
    return logger