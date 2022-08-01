import telebot
from requests import get
   
bot = telebot.TeleBot("******************************************************")

@bot.message_handler(commands = ['start'])
def start(message):
    mess = f'Привет, <b>{message.from_user.first_name} <u>{message.from_user.last_name}</u></b>'
    bot.send_message(message.chat.id, mess, parse_mode='html')

@bot.message_handler(commands = ['bazis'])
def bazis(message):
    bazis = open('Презентация Базис.pdf', 'rb')
    bot.send_document(message.chat.id, bazis)

@bot.message_handler(commands = ['finrez'])
def finrez(message):
    finrez = open('Презентация Фин резерв.pdf', 'rb')
    bot.send_document(message.chat.id, finrez)

@bot.message_handler(commands = ['rules'])
def rules(message):
    rules = open('rules.pdf', 'rb')
    bot.send_document(message.chat.id, rules)

@bot.message_handler(commands = ['license'])
def license(message):
    license = open('Лицензия.pdf', 'rb')
    bot.send_document(message.chat.id, license)

@bot.message_handler(commands = ['zdozitie'])
def zdozitie(message):
    zdozitie = open('Заявление на дожитие.pdf', 'rb')
    bot.send_document(message.chat.id, zdozitie)

@bot.message_handler(commands = ['znefin'])
def znefin(message):
    znefin = open('Заявление нефин изменения.pdf', 'rb')
    bot.send_document(message.chat.id, znefin)

@bot.message_handler(commands = ['zstrah'])
def zstrah(message):
    zstrah = open('Заявление на страхование.pdf', 'rb')
    bot.send_document(message.chat.id, zstrah)

@bot.message_handler(commands = ['zoff'])
def zoff(message):
    zoff = open('Заявление на расторжение.pdf', 'rb')
    bot.send_document(message.chat.id, zoff)

@bot.message_handler(commands = ['photo'])
def photo(message):
#    photo = open('photo.png', 'rb')
    bot.send_photo(message.chat.id, get('http://thecatapi.com/api/images/get?format=src').content)

@bot.message_handler()
def get_user_text(message):
    if message.text == "Hello":
        bot.send_message(message.chat.id, "And you", parse_mode='html')
    elif message.text == "id":
        bot.send_message(message.chat.id, f"Your ID {message.from_user.id}", parse_mode='html')
    elif message.text == "photo":
        photo = open('photo.png', 'rb')
        bot.send_photo(message.chat.id, photo)
    else:
        bot.send_message(message.chat.id, "Выберите пункт меню", parse_mode='html')

bot.polling(none_stop=True)