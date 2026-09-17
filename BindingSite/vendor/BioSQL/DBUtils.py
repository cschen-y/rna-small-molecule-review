











"""Helper code for Biopython's BioSQL code (for internal use)."""

import os

_dbutils: dict[str, type["Generic_dbutils"]] = {}


class Generic_dbutils:
    """Default database utilities."""

    def __init__(self):
        """Create a Generic_dbutils object."""

    def tname(self, table):
        """Return the name of the table."""
        if table != "biosequence":
            return table
        else:
            return "bioentry"

    def last_id(self, cursor, table):
        """Return the last used id for a table."""
        
        table = self.tname(table)
        sql = f"select max({table}_id) from {table}"
        cursor.execute(sql)
        rv = cursor.fetchone()
        return rv[0]

    def execute(self, cursor, sql, args=None):
        """Just execute an sql command."""
        cursor.execute(sql, args or ())

    def executemany(self, cursor, sql, seq):
        """Execute many sql commands."""
        cursor.executemany(sql, seq)

    def autocommit(self, conn, y=1):
        """Set autocommit on the database connection."""
        


class Sqlite_dbutils(Generic_dbutils):
    """Custom database utilities for SQLite."""

    def _sub_placeholder(self, sql):
        """Format the argument placeholders for sqlite (PRIVATE)."""
        return sql.replace("%s", "?")

    def execute(self, cursor, sql, args=None):
        """Execute SQL command.

        Replaces %s with ? for variable substitution in sqlite3.
        """
        sql = self._sub_placeholder(sql)
        cursor.execute(sql, args or ())

    def executemany(self, cursor, sql, seq):
        """Execute many sql statements."""
        sql = self._sub_placeholder(sql)
        cursor.executemany(sql, seq)


_dbutils["sqlite3"] = Sqlite_dbutils


class Mysql_dbutils(Generic_dbutils):
    """Custom database utilities for MySQL."""

    def last_id(self, cursor, table):
        """Return the last used id for a table."""
        if os.name == "java":
            return Generic_dbutils.last_id(self, cursor, table)
        try:
            
            return cursor.insert_id()
        except AttributeError:
            
            
            
            return cursor.lastrowid


_dbutils["MySQLdb"] = Mysql_dbutils


class _PostgreSQL_dbutils(Generic_dbutils):
    """Base class for any PostgreSQL adaptor."""

    def next_id(self, cursor, table):
        table = self.tname(table)
        sql = f"SELECT nextval('{table}_pk_seq')"
        cursor.execute(sql)
        rv = cursor.fetchone()
        return rv[0]

    def last_id(self, cursor, table):
        table = self.tname(table)
        sql = f"SELECT currval('{table}_pk_seq')"
        cursor.execute(sql)
        rv = cursor.fetchone()
        return rv[0]


class Psycopg2_dbutils(_PostgreSQL_dbutils):
    """Custom database utilities for Psycopg2 (PostgreSQL)."""

    def autocommit(self, conn, y=True):
        """Set autocommit on the database connection."""
        if y:
            if os.name == "java":
                conn.autocommit = 1
            else:
                conn.set_isolation_level(0)
        else:
            if os.name == "java":
                conn.autocommit = 0
            else:
                conn.set_isolation_level(1)


_dbutils["psycopg2"] = Psycopg2_dbutils


class Pgdb_dbutils(_PostgreSQL_dbutils):
    """Custom database utilities for Pgdb (aka PyGreSQL, for PostgreSQL)."""

    def autocommit(self, conn, y=True):
        """Set autocommit on the database connection. Currently not implemented."""
        raise NotImplementedError("pgdb does not support this!")


_dbutils["pgdb"] = Pgdb_dbutils


def get_dbutils(module_name):
    """Return the correct dbutils object for the database driver."""
    try:
        return _dbutils[module_name]()
    except KeyError:
        return Generic_dbutils()
