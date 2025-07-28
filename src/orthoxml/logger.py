# logger.py

import logging

# Configure the base logger
logger = logging.getLogger("orthoxml")
logging_level = logging.WARNING

# Create formatters
console_formatter = logging.Formatter(
    fmt='%(asctime)s - %(module)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# Console handler
console_handler = logging.StreamHandler()
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
    base = "orthoxml"
    full_name = f"{base}.{name}" if name else base
    return logging.getLogger(full_name)

def set_logger_level(level):
    level = getattr(logging, level.upper(), logging.WARNING)

    logger.setLevel(level)
    for handler in logger.handlers:
        handler.setLevel(level)