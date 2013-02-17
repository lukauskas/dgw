import io


class ChunkedIO(io.IOBase):
    """
    Wrapper around file object that would chunk read and write operations
    as to fix any nonsenses that tend to happen due to incorrect handling
    of large objects.

    see Issues #2806, #574 for instance
    """
    __max_chunk_size = None
    __file = None

    def __init__(self, file, max_chunk_size=2 ** 28):
        """
        Initialise `ChunkedIO` object with a file_like buffer.

        :param file:
        :param max_chunk_size: the maximum size of the chunks that will be read from or written to the file provided
        :return:
        """
        max_chunk_size = int(max_chunk_size)
        if max_chunk_size <= 0:
            raise ValueError("Max chunk size should be greater than 0")

        self.__max_chunk_size = max_chunk_size
        self.__file = file

    def read(self, size=None):
        max_chunk_size = self.__max_chunk_size
        if size is not None and size < max_chunk_size:
            return self.__file.read(size)

        data_read = ""
        while size is None or size > 0:
            read_size = max_chunk_size if size is None or size >= max_chunk_size else size
            current_read = self.__file.read(size)
            data_read += current_read
            # If EOF was reached
            if len(current_read) < read_size:
                break
            elif size is not None:
                size -= len(current_read)

        return data_read

    def write(self, p_str):
        max_chunk_size = self.__max_chunk_size

        while len(p_str) > 0:
            current_chunk_size = min(max_chunk_size, len(p_str))
            current_chunk = p_str[:current_chunk_size]
            p_str = p_str[current_chunk_size:]

            self.__file.write(current_chunk)

    def seek(self, *args, **kwargs):
        self.file.seek(*args, **kwargs)

    def close(self, *args, **kwargs):
        self.__file.close(*args, **kwargs)

