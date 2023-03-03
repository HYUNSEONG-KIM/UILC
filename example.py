from uilc import 



s = input("Enter Lambertian number: ")
if s.isnumeric():
    s = float(s)
if s <0:
    raise ValueError("s must be larger than 0")

