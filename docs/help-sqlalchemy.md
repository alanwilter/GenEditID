# SQLAlchemy

## Introduction
SQLAlchemy is a deep and powerful thing made up of many layers. This cheat sheet sticks to parts of the ORM (Object Relational Mapper) layer, and aims to be a reference not a tutorial. That said, if you are familiar with SQL then this cheat sheet should get you well on your way to understanding SQLAlchemy.

## Basic Models

One model is used to describe one database table. For example:

```Python
from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy import (
    Column,
    Integer,
    String,
    Boolean,
    ForeignKey,
    DateTime,
    Sequence,
    Float
)
import datetime

Base = declarative_base()


class Book(Base):  #<-------------------------
    __tablename__  = "books"    # matches the name of the actual database table
    id             = Column(Integer, Sequence('book_seq'), primary_key=True) # plays nice with all major database engines
    name           = Column(String(50))                                      # string column need lengths
    author_id      = Column(Integer, ForeignKey('authors.id'))               # assumes there is a table in the database called 'authors' that has an 'id' column
    price          = Column(Float)
    date_added     = Column(DateTime, default=datetime.datetime.now)         # defaults can be specified as functions
    promote        = Column(Boolean,default=False)                           #     or as values

```

## Queries and Interactions

```Python
from sqlalchemy import create_engine
engine = create_engine(DATABASE_URI)

from sqlalchemy.orm import sessionmaker
DBsession = sessionmaker(bind=engine)
DBsession = session()

# fetch everything
lBooks = DBSession.query(Book)  # returns a Query object.
for oBook in lBooks:
    print oBook.name

# simple filters
lBooks = DBSession.query(Book).filter_by(author_id=1) #returns all the books for a specific author

# more complex filters
lBooks = DBSession.query(Book).filter(Book.price<20) #returns all the books with price <20. Note we use filter, not filter_by

# filters can be combined
lBooks = DBSession.query(Book).filter_by(author_id=1).filter(Book.price<20) #all books by a specific author, with price<20

# logical operations can be used in filters
from sqlalchemy import or_
lBooks = DBSession.query(Book).filter(or_(Book.price<20,promote==True)) # returns all books  that cost less than 20 OR are being promoted

# ordering
from sqlalchemy import desc
DBSession.query(Book).order_by(Book.price) #get all books ordered by price
DBSession.query(Book).order_by(desc(Book.price)) #get all books ordered by price descending

# other useful things
DBSession.query(Book).count() #returns the number of books
DBSession.query(Book).offset(5) #offset the result by 5
DBSession.query(Book).limit(5) # return at most 5 books
DBSession.query(Book).first() #return the first book only or None
DBSession.query(Book).get(8) #return the Book with primary key = 8, or None
```

## Relationships

Relationships between SQL tables are described in terms of foreign key relationships. From the example, the `books` table has a foreign key field pointing to the `id` field of the `authors` table. SQLAlchemy makes leveraging and examining those relationships pretty straight forward.

### One to many relationships

Assume we are keeping track of the books of various authors. One author can have many books.

```Python
class Book(Base):
    __tablename__  = "books"    # matches the name of the actual database table
    id             = Column(Integer, Sequence('book_seq'), primary_key=True)
    name           = Column(String(50))                                    
    author_id      = Column(Integer, ForeignKey('authors.id'))              
    author         = relationship("Author", backref="books")  # <-----------------------

class Author(Base):
    __tablename__  = "books"    # matches the name of the actual database table
    id             = Column(Integer, Sequence('book_seq'), primary_key=True)
    name           = Column(String(50))
```

The marked line configures the relationship between the models. Note that "Author" is a string. It doesn't need to be, it can also be a class. Using a string here removes the possibility of certain NameErrors. Note that the relationship is configured in both directions in one line. A book's author is accessible via the author attribute, and an author's books are accessible via the author's books attribute.

Here are a few ways you can make use of the relationship once it is configured:

```Python
oBook = DBSession.query(Book).filter_by(name="Harry Potter and the methods of rationality").first()
oAuthor = oBook.author   # oAuthor is now an Author instance. oAuthor.id == oBook.author_id

# it works the other way as well
oAuthor = DBSession.query(Author).filter_by(name="Orsan Scott Card")
for oBook in oAuthor.books:
    print oBook.name

# adding a new book
oNewBook = Book()
oBook.name = "Ender's Game"
oBook.author = oAuthor

# adding a new book in a different way...
oNewBook = Book()
oBook.name = "Ender's Shadow"
oAuthor.books.append(oBook)
```

## Engine configuration

### Connection Strings

```Python
#t he general form of a connection string:
`dialect+driver://username:password@host:port/database`

# SQLITE:
'sqlite:///:memory:' # store everything in memory, data is lost when program exits
'sqlite:////absolute/path/to/project.db')  # Unix/Mac
'sqlite:///C:\\path\\to\\project.db' # Windows
r'sqlite:///C:\path\to\project.db' # Windows alternative

# PostgreSQL
'postgresql://user:pass@localhost/mydatabase'
'postgresql+psycopg2://user:pass@localhost/mydatabase'
'postgresql+pg8000://user:pass@localhost/mydatabase'
```

### Engine, Session and Base

```Python
# set up the engine
engine = create_engine(sConnectionString, echo=True)   
# echo=True makes the sql commands issued by sqlalchemy get output to the console, useful for debugging

# bind the dbsession to the engine
DBSession.configure(bind=engine)

# now you can interact with the database if it exists

# import all your models then execute this to create any tables that don't yet exist. This does not handle migrations
Base.metadata.create_all(engine)       
```

## Reference

[A SQLAlchemy Cheat Sheet](https://www.codementor.io/sheena/understanding-sqlalchemy-cheat-sheet-du107lawl)
