#!/usr/bin/env python

from __future__ import print_function
import pika
from time import sleep, time
import timeout_decorator

@timeout_decorator.timeout(5)
def rmq_listen():
    channel.start_consuming()

# Part 1: send message
print("Start sending message...")
connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
channel = connection.channel()
channel.queue_declare(queue='hello')
msg = 'Hello World! {}'.format(time())
channel.basic_publish(exchange='',routing_key='hello', body=msg)
print(" [x] Sent '{}'".format(msg))
connection.close()

print("done sending message. Wait for 5 seconds and start receiving...")
sleep(5)

# Part 2: receive message
connection = pika.BlockingConnection(pika.ConnectionParameters(host='localhost'))
channel = connection.channel()
channel.queue_declare(queue='hello')
print(' [*] Waiting 5 seconds for messages.')
def callback(ch, method, properties, body):
    print(" [x] Received %r" % (body,))
channel.basic_consume(queue='hello',on_message_callback=callback,auto_ack=False)
try:
    rmq_listen()
except timeout_decorator.timeout_decorator.TimeoutError:
    print("Done waiting, now close connection")
    connection.close()
