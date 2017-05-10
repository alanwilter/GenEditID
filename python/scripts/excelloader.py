from datetime import datetime


class ExcelLoader:

    def get_value(self, value):
        if not value:
            return None
        elif str(value) == 'nan':
            return None
        elif value == '':
            return None
        return value

    def get_string(self, value, max_length=-1):
        if self.get_value(value):
            value = str(value).strip()
            if max_length >= 0 and len(value) > max_length:
                if max_length > 20:
                    return value[:max_length - 3] + "..."
                else:
                    return value[:max_length]
            return value

    def get_int(self, value):
        if self.get_value(value):
            return int(value)

    def get_float(self, value):
        if self.get_value(value):
            return float(value)

    def get_date(self, value, format='%Y%m%d'):
        if self.get_value(value):
            return datetime.strptime(str(value), '%Y%m%d')
