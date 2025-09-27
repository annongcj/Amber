#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#define _TEST_HASH_STANDALONE 1


/* Fast and simple hash table 
   used for input processing.
*/
#include "hash.h"


unsigned long hash_function(char* str) {
    unsigned long i = 0;
    int           j;
    for (j=0; ((char *)str)[j]; j++) {
        i  = (( i + str[j] ) * 11) % CAPACITY;
    }
    return( i );
}
 
 
static LinkedList* allocate_list () {
    // Allocates memory for a Linkedlist pointer
    LinkedList* list = (LinkedList*) malloc (sizeof(LinkedList));
    return( list );
}
 
static LinkedList* linkedlist_insert(LinkedList* list, Ht_item* item) {
    // Inserts the item onto the Linked List
    if (!list) {
        LinkedList* head = allocate_list();
        head->item = item;
        head->next = NULL;
        list = head;
        return list;
    } 
     
    else if (list->next == NULL) {
        LinkedList* node = allocate_list();
        node->item = item;
        node->next = NULL;
        list->next = node;
        return list;
    }
 
    LinkedList* temp = list;
    while (temp->next->next) {
        temp = temp->next;
    }
     
    LinkedList* node = allocate_list();
    node->item = item;
    node->next = NULL;
    temp->next = node;
     
    return list;
}
 
static Ht_item* linkedlist_remove(LinkedList* list) {
    // Removes the head from the linked list
    // and returns the item of the popped element
    if (!list)
        return NULL;
    if (!list->next)
        return NULL;
    LinkedList* node = list->next;
    LinkedList* temp = list;
    temp->next = NULL;
    list = node;
    Ht_item* it = NULL;
    memcpy(temp->item, it, sizeof(Ht_item));
    free(temp->item->key);
    free(temp->item->value);
    free(temp->item);
    free(temp);
    return it;
}
 
static void free_linkedlist(LinkedList* list) {
    LinkedList* temp = list;
    while (list) {
        temp = list;
        list = list->next;
        free(temp->item->key);
        free(temp->item->value);
        free(temp->item);
        free(temp);
    }
}
 
static LinkedList** create_overflow_buckets(HashTable* table) {
    int i;

    // Create the overflow buckets; an array of linkedlists
    LinkedList** buckets = (LinkedList**) calloc (table->size, sizeof(LinkedList*));
    for (i=0; i<table->size; i++)
        buckets[i] = NULL;
    return buckets;
}
 
static void free_overflow_buckets(HashTable* table) {
    int i;

    // Free all the overflow bucket lists
    LinkedList** buckets = table->overflow_buckets;
    for (i=0; i<table->size; i++)
        free_linkedlist(buckets[i]);
    free(buckets);
}
 
 
Ht_item* create_item(char* key, char* value) {
    // Creates a pointer to a new hash table item
    Ht_item* item = (Ht_item*) malloc (sizeof(Ht_item));
    item->key = (char*) malloc ( sizeof(char) * (strlen(key) + 1));
    item->value = (char*) malloc (sizeof(char) * (strlen(value) + 1));
     
    strcpy(item->key, key);
    strcpy(item->value, value);
 
    return item;
}
 
HashTable* create_table(int size) {
    int i;

    // Creates a new HashTable
    HashTable* table = (HashTable*) malloc (sizeof(HashTable));
    
    table->size = size;
    table->count = 0;
    table->items = (Ht_item**) calloc (table->size, sizeof(Ht_item*));
    for (i=0; i<table->size; i++)
        table->items[i] = NULL;
    table->overflow_buckets = create_overflow_buckets(table);
    table->keyList          = NULL;
 
    return table;
}
 
void free_item(Ht_item* item) {
    // Frees an item
    free(item->key);
    free(item->value);
    free(item);
}
 
void free_table(HashTable* table) {
    // Frees the table
    int i;
    for (i=0; i<table->size; i++) {
        Ht_item* item = table->items[i];
        if (item != NULL)
            free_item(item);
    }
 
    free_overflow_buckets(table);
    free(table->items);
    free(table);
}
 
void handle_collision(HashTable* table, unsigned long index, Ht_item* item) {
    LinkedList* head = table->overflow_buckets[index];
 
    if (head == NULL) {
        // We need to create the list
        head = allocate_list();
        head->item = item;
        table->overflow_buckets[index] = head;
        return;
    }
    else {
        // Insert to the list
        table->overflow_buckets[index] = linkedlist_insert(head, item);
        return;
    }
 }
 
void ht_insert(HashTable* table, char* key, char* value) {
    // Create the item
    Ht_item    *item = create_item(key, value);
    LinkedList *keyForList;
 
    // Compute the index
    unsigned long index = hash_function(key);
 
    Ht_item* current_item = table->items[index];
     
    if (current_item == NULL) {
        // Key does not exist.
        if (table->count == table->size) {
            // Hash Table Full
            printf("Insert Error: Hash Table is full\n");
            // Remove the failed item
            free_item(item);
            return;
        }
         
        // Insert directly
        table->items[index] = item; 
        table->count++;
        
        // Keep a list of items that have been added to the table:
        keyForList = (LinkedList *)allocate_list();
        keyForList->item     = item;  //just point to the hash item
        keyForList->next     = table->keyList; //insert at head of list.
        table->keyList       = keyForList; 
    }
    else {
        // Overwrite: We only need to update value
        if (strcmp(current_item->key, key) == 0) {
                printf("Hash insert Warning: overwriting entry .%s.\n", key);
                strcpy(table->items[index]->value, value);
                return;
        }
        // Collision
        else {
            printf("Hash insert Warning: collision .%s.\n", key);
            handle_collision(table, index, item);

            // Keep a list of items that have been added to the table:
            keyForList = (LinkedList *)allocate_list();
            keyForList->item     = item;  //just point to the hash item
            keyForList->next     = table->keyList; //insert at head of list.
            table->keyList       = keyForList; 
            return;
        }
    }
}
 
char* ht_search(HashTable* table, char* key) {
    // Searches the key in the hashtable
    // and returns NULL if it doesn't exist
    int index = hash_function(key);
    Ht_item* item = table->items[index];
    LinkedList* head = table->overflow_buckets[index];
 
    // Ensure that we move to items which are not NULL
    while (item != NULL) {
        if (strcmp(item->key, key) == 0)
            return item->value;
        if (head == NULL)
            return NULL;
        item = head->item;
        head = head->next;
    }
    return NULL;
}
 
void ht_delete(HashTable* table, char* key) {
    // Deletes an item from the table
    int        index = hash_function(key);
    Ht_item*    item = table->items[index];
    LinkedList* head = table->overflow_buckets[index];
 
    if (item == NULL) {
        // Does not exist. Return
        return;
    }
    else {
        if (head == NULL && strcmp(item->key, key) == 0) {
            // No collision chain. Remove the item
            // and set table index to NULL

            table->items[index] = NULL;

            //clear from the list of members, linear in time.
            LinkedList *keyList_now, *keyList_prev;            
            keyList_prev = NULL;
            keyList_now  = table->keyList;
            while( keyList_now ){
                //snip out this item?
                if( strcmp(keyList_now->item->key, key) == 0 ){
                    if( keyList_now == table->keyList ){ 
                         table->keyList = keyList_now->next;
                    }else{
                         keyList_prev->next  = keyList_now->next;
                    } 
                    break;
                }
                keyList_prev = keyList_now;
                keyList_now  = keyList_now->next;
            } 

            free_item(item);
            table->count--;
            return;
        }
        else if (head != NULL) {
            // Collision Chain exists
            if (strcmp(item->key, key) == 0) {
                //clear from the list of members, linear in time.
                LinkedList *keyList_now, *keyList_prev;            
                keyList_prev = NULL;
                keyList_now  = table->keyList;
                while( keyList_now ){
                    //snip out this item?
                    if( strcmp(keyList_now->item->key, key) == 0 ){
                        if( keyList_now == table->keyList ){ 
                             table->keyList = keyList_now->next;
                        }else{
                             keyList_prev->next  = keyList_now->next;
                        } 
                        break;
                    }
                    keyList_prev = keyList_now;
                    keyList_now  = keyList_now->next;
                } 


                // Remove this item and set the head of the list
                // as the new item
                 
                free_item(item);
                LinkedList* node = head;
                head = head->next;
                node->next = NULL;
                table->items[index] = create_item(node->item->key, node->item->value);
                free_linkedlist(node);
                table->overflow_buckets[index] = head;
                return;
            }
 
            LinkedList* curr = head;
            LinkedList* prev = NULL;
             
            while (curr) {
                if (strcmp(curr->item->key, key) == 0) {
                    if (prev == NULL) {
                        // First element of the chain. Remove the chain
                        free_linkedlist(head);
                        table->overflow_buckets[index] = NULL;
                        return;
                    }
                    else {
                        // This is somewhere in the chain
                        prev->next = curr->next;
                        curr->next = NULL;
                        free_linkedlist(curr);
                        table->overflow_buckets[index] = head;
                        return;
                    }
                }
                curr = curr->next;
                prev = curr;
            }
        }
    }
}
 
void print_search(HashTable* table, char* key) {
    char* val;
    if ((val = ht_search(table, key)) == NULL) {
        printf("%s does not exist\n", key);
        return;
    }
    else {
        printf("Key:%s, Value:%s\n", key, val);
    }
}
 
void print_table(HashTable* table) {

    int i;

    printf("\n-------------------\n");
    for (i=0; i<table->size; i++) {
        if (table->items[i]) {
            printf("Index:%d, Key:%s, Value:%s", i, table->items[i]->key, table->items[i]->value);
            if (table->overflow_buckets[i]) {
                printf(" => Overflow Bucket => ");
                LinkedList* head = table->overflow_buckets[i];
                while (head) {
                    printf("Key:%s, Value:%s ", head->item->key, head->item->value);
                    head = head->next;
                }
            }
            printf("\n");
        }
    }
    printf("-------------------\n");
}

#ifdef TEST_HASH_STANDALONE

int main() {
    HashTable* ht = create_table(CAPACITY);
    ht_insert(ht, "11", "First address");
    ht_insert(ht, "22", "Second address");
    ht_insert(ht, "Sausage", "Third address");
    ht_insert(ht, "Beenz",   "Fourth address");
    print_search(ht, "11");
    print_search(ht, "22");
    print_search(ht, "3");
    print_search(ht, "Sausage");
    print_search(ht, "Beenz");  
    print_table(ht);
    ht_delete(ht, "11");
    ht_delete(ht, "22");
    print_table(ht);
    free_table(ht);
    return 0;
}

#endif

