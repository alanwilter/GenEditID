from datetime import datetime
import pandas


class ExcelLoader:

    def get_value(self, text):
        if not text:
            return None
        elif pandas.isnull(text):
            return None
        elif str(text) == 'nan':
            return None
        elif text == '':
            return None
        return text

    def get_string(self, text, max_length=-1):
        if self.get_value(text):
            value = str(text).strip()
            if max_length >= 0 and len(value) > max_length:
                if max_length > 20:
                    return value[:max_length - 3] + "..."
                else:
                    return value[:max_length]
            return value

    def get_int(self, text):
        if self.get_value(text):
            try:
                return int(text)
            except ValueError:
                return None
        return None

    def get_float(self, text):
        if self.get_value(text):
            try:
                return float(text)
            except ValueError:
                return None
        return None

    def get_date(self, text, format='%Y%m%d'):
        if self.get_value(text):
            try:
                return datetime.strptime(str(text), '%Y%m%d')
            except ValueError:
                return None
        return None
