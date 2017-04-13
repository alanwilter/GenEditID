from datetime import datetime


class ExcelLoader:

    def get_string(self, cell, max_length=-1):
        if not cell:
            return None

        if cell.value is None:
            return None

        value = str(cell.value).strip()

        if value == '':
            return None

        if max_length >= 0 and len(value) > max_length:
            if max_length > 20:
                return value[:max_length - 3] + "..."
            else:
                return value[:max_length]

        return value

    def get_int(self, cell):
        if not cell:
            return None

        if cell.value is None:
            return None

        if isinstance(cell.value, int):
            return cell.value

        return int(str(cell.value))

    def get_float(self, cell):
        if not cell:
            return None

        if cell.value is None:
            return None

        if isinstance(cell.value, float):
            return cell.value

        return float(str(cell.value))

    def get_date(self, cell, format='%Y%m%d'):
        if not cell:
            return None

        if cell.value is None:
            return None

        if isinstance(cell.value, datetime):
            return cell.value

        return datetime.strptime(str(cell.value), '%Y%m%d')
