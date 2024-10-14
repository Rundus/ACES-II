

x = 5

def func():
    global x
    x += 5

for i in range(5):
    func()
    print(x)
