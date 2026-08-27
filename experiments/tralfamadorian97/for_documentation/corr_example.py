import numpy as np



def go():
    example_1 = np.array([[1,1,1],
                         [1,1,0],
                         [1,0,1],
                          [1,0,0],
                          [0,1,1],
                          [0,1,0],
                          [0,0,1],
                          [0,0,0]
                         ]
                         )

    result=np.corrcoef(example_1, rowvar=False)
    example_2 = np.array([
                          [1,1,0],
                          [1,0,1],
                          [0,1,1],
                          [0,0,0]
                          ]

                         )
    result_2=np.corrcoef(example_2, rowvar=False)
    print(result)
    print(result_2)

if __name__ == '__main__':
    go()