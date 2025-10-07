import logs


class EmptyParsingError(Exception):
    # If nothing is extracted from an external file
    def __init__(self, msg):
        # Call the base class constructor with the parameters it needs
        super().__init__(msg)
        logs.write_in_log(msg, "ERROR")

    pass


class DimError(Exception):
    # If dimensions of matrix and vectors are not consistent
    def __init__(self, msg):
        # Call the base class constructor with the parameters it needs
        super().__init__(msg)
        logs.write_in_log(msg, "ERROR")

    pass


class SensInputError(Exception):
    def __init__(self, msg):
        # Call the base class constructor with the parameters it needs
        super().__init__(msg)
        logs.write_in_log(msg, "ERROR")

    pass


class UserInputError(Exception):
    # If user input data is of the wrong type
    def __init__(self, msg):
        # Call the base class constructor with the parameters it needs
        super().__init__(msg)
        logs.write_in_log(msg, "ERROR")

    pass


class MissingDataError(Exception):
    # If requested data doesnt exist
    def __init__(self, msg):
        # Call the base class constructor with the parameters it needs
        super().__init__(msg)
        logs.write_in_log(msg, "ERROR")

    pass
