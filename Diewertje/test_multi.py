import numpy as np
from time import time , sleep
import multiprocessing as mp # Here we add the multiprocessing package !


### We leave thi s part of the program untouched ###
###################################################
# Function that takes a number and returns e ^[number ] .
def exponent ( a ) :
    ans = np. exp( a ) # Calculate exponent # Wait f i v e seconds
    return ans # Return exponent

# Function that takes a number and returns 1 / [number ] .
def divide (b) :
    ans = 1.0/b # Calculate exponent # Wait f i v e seconds
    return ans # Return exponent



def worker_func(inQ,outQ) :

    # Inf ini t e loop to keep the function running
    while True :
        # t r y / except cons truct ion to catch any e r ror s .
        try :
            
            input = inQ.get()
            
            if input == None:
                break

# Spl i t the input into a number and a function
            number = input[0]
            function = input[1]
# Calculate the r e sul t
            output = [function(number),function]

            outQ.put(output)
        
        except (Exception) as e :
            print( 'Error!' , e.message)
            break
        
def main():  
    # Jus t some values to cal culat e the exponents of
    a = np.log( 3.1 )
    b = 1.0/2.4
     
    inputQ = mp.Queue( )
    outputQ = mp.Queue( )
    
    Workers = [mp.Process(target=worker_func, args =(inputQ, outputQ)),
               mp.Process(target=worker_func, args =(inputQ, outputQ))]
    
    Workers[0].start()
    Workers[1].start()
    
    
    input1 = [a, exponent]
    input2 = [b, divide]
    
    time_one = time()
    
    inputQ.put(input1)
    inputQ.put(input2)
    
    result1 = outputQ.get()
    result2 = outputQ.get()
    
    # Stop the timer
    time_two = time ( )
    
    # Close the worker proces ses
    # Fi r s t stop the function from running
    inputQ.put(None)
    inputQ.put(None)
    
    # Then c los e and join the proces ses
    for w in Workers :
        w.join()
    # Print the r e sul t s
    print(result1, result2)
    
if __name__ == "__main__": 
    main()