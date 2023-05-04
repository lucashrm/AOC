#include <stdlib.h>
#include <semaphores.h>

int main(void) { 
    sem_timedwait(NULL, NULL);
    return 0; 
}

