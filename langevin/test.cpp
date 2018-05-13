#include<iostream>
#include<vector>

int main(){
  std::vector<int> myvector;
  int myint;
  char str[10];
  printf("Please enter some integers (enter 0 to end):\n");
  do {
    scanf("%d",&myint);
    myvector.push_back(myint);
  } while (myint);
  printf("myvector stores %d numbers.\n",int(myvector.size()));
  for(int i=0;i<int(myvector.size());i++){
    printf("%d\n",myvector[i]);
  }

  return 0;
}
    
