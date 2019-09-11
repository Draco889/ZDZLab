import mysql.connector

mydb = mysql.connector.connect(
  host="127.0.0.1",
  user="Mark",
  passwd="Draco889Tiger131",
  database="snp"
)

mycursor = mydb.cursor()

mycursor.execute("DROP TABLE position")

TABLE = (
    "CREATE TABLE `position` ("
    "  `rsid` varchar(50) NOT NULL,"
    "  `chromosome` int(2) NOT NULL,"
    "  `bp` int(25) NOT NULL,"
    "  PRIMARY KEY (`rsid`)"
    ") ENGINE=InnoDB")

mycursor.execute(TABLE)

PATH = "/Users/Mark/PycharmProjects/ZDZLab1/Reference/"

o = open(PATH + "homo_sapiens-chr1.vcf", "r")

add_record = ("INSERT INTO position "
               "(rsid, chromosome, bp) "
               "VALUES (%s, %s, %s)")

for line in o:
    if "#" in line:
        continue
    else:
        line = line.split()
        record = (line[2], int(line[0]), int(line[1]))
        mycursor.execute(add_record, record)

mydb.commit()
mycursor.close()
mydb.close()

