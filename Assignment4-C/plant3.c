#include "include/pbPlots.c"
#include "include/supportLib.c"
#include "include/arrayString.c"
#include "include/arrayFloat.c"
#include "include/array.c"

#define TIME 2000
#define EPOCH 500
#define LEARNING_RATE_END 0.01

// parameter jizz: EPOCH 500, LEARNING RATE 0.1, LEARNING RATE END 0.01 
// parameter tanpa epoch: EPOCH 1, LEARNING RATE 0.16
float LEARNING_RATE = 0.1;
typedef struct {
    ArrayString membership;
    ArrayFloat miu;
} fuzzifier;

float sgn(float val) {
  if (val > 0.0) {
    return 1.0;
  } else if (val < 0.0) {
    return -1.0;
  } else {
    return 0;
  }
}

ArrayString rule(ArrayString* membership_ek, ArrayString* delta_membership_ek, ArrayFloat* rule){
    ArrayString rules;
    initArrayString(&rules, 1);

    for (size_t i = 0; i < membership_ek->used; i++){
        for (size_t j = 0; j < delta_membership_ek->used; j++){
            if((strcmp(membership_ek->array[i], "N") == 0) && (strcmp(delta_membership_ek->array[j], "N") == 0)){
                insertArrayString(&rules, "N");
                rule->array[0] = 1.0;
                // printf("i: %d, j: %d", i, j);
            }
            if((strcmp(membership_ek->array[i], "N") == 0) && (strcmp(delta_membership_ek->array[j], "Z") == 0)){
                insertArrayString(&rules, "N");
                rule->array[3] = 1.0;
                // printf("i: %d, j: %d", i, j);
            }
            if((strcmp(membership_ek->array[i], "N") == 0) && (strcmp(delta_membership_ek->array[j], "P") == 0)) {
                insertArrayString(&rules, "N");
                rule->array[6] = 1.0;
                // printf("i: %d, j: %d", i, j);
            }
                
            
            if((strcmp(membership_ek->array[i], "Z") == 0) && (strcmp(delta_membership_ek->array[j], "N") == 0)){
                insertArrayString(&rules, "N");
                rule->array[1] = 1.0;
                // printf("i: %d, j: %d", i, j);
            }
            if((strcmp(membership_ek->array[i], "Z") == 0) && (strcmp(delta_membership_ek->array[j], "Z") == 0)){
                insertArrayString(&rules, "Z");
                rule->array[4] = 1.0;
                // printf("i: %d, j: %d", i, j);
            }
            if((strcmp(membership_ek->array[i], "Z") == 0) && (strcmp(delta_membership_ek->array[j], "P") == 0)) {
                insertArrayString(&rules, "P");
                rule->array[7] = 1.0;
                // printf("i: %d, j: %d", i, j);
            }
            
            if((strcmp(membership_ek->array[i], "P") == 0) && (strcmp(delta_membership_ek->array[j], "N") == 0)){
                insertArrayString(&rules, "P");
                rule->array[2] = 1.0;
                // printf("i: %d, j: %d", i, j);
            }
            if((strcmp(membership_ek->array[i], "P") == 0) && (strcmp(delta_membership_ek->array[j], "Z") == 0)){
                insertArrayString(&rules, "P");
                rule->array[5] = 1.0;
                // printf("i: %d, j: %d", i, j);
            }
            if((strcmp(membership_ek->array[i], "P") == 0) && (strcmp(delta_membership_ek->array[j], "P") == 0)) {
                insertArrayString(&rules, "P");
                rule->array[8] = 1.0;
                // printf("i: %d, j: %d", i, j);
            }
            
        }
    }
    return rules;
}

fuzzifier fuzzifierFunc(float ek, float a, float b, float c){
    fuzzifier data;

    initArrayString(&data.membership, 1);
    initArrayFloat(&data.miu, 1);

    float miu_temp = 0.0;

    float a_ref = a;
    float b_ref = b;
    float c_ref = c;

    // gradient
    float m_N = (0-1)/(b_ref - a_ref);

    float m_Z1 = (1-0)/(b_ref - a_ref);
    float m_Z2 = (0-1)/(c_ref - b_ref);

    float m_P = (1-0)/(c_ref - b_ref);

    // Calculating MIU
    if(ek<= a_ref || (ek>a_ref && ek<b_ref)){
        insertArrayString(&data.membership, "N");
        if(ek<= a_ref){
            miu_temp = 1.0;
        }
        else if(ek>a_ref && ek<b_ref){
            miu_temp = m_N*(ek -a_ref) + 1;
        }
        insertArrayFloat(&data.miu, miu_temp);
    }
        
    if(ek> a_ref && ek<c_ref){
        insertArrayString(&data.membership, "Z");
        if(ek> a_ref && ek<b_ref){
            miu_temp = m_Z1*(ek -a_ref) + 0;
        }
        else{
            miu_temp = m_Z2*(ek -b_ref) + 1;
        }
        // printf("kon: %f\n", miu_temp);
        insertArrayFloat(&data.miu, miu_temp);
    }
    if(ek> b_ref){
        insertArrayString(&data.membership, "P");
        if(ek>=c_ref){
            miu_temp = 1.0;
        }
        else{
            miu_temp = m_P*(ek-b_ref) + 0;
        }
        insertArrayFloat(&data.miu, miu_temp);
    }

    return data;
    
}

ArrayFloat defuzzification(ArrayFloat *miu_ek, ArrayFloat *miu_delta_ek, ArrayFloat *rules){
    ArrayFloat weight;
    initArrayFloat(&weight, 1);
    float weight_temp = 0.0;

    int index_e=0, index_delta_e=0;

    for(int m=0;m<9;m++){
        insertArrayFloat(&weight, 0);
    }

    for(int m=0;m<9;m++){
        if(index_delta_e < miu_delta_ek->used){
            if(rules->array[m] > 0){
                weight_temp = fmin(miu_ek->array[index_e], miu_delta_ek->array[index_delta_e]);
                weight.array[m] = weight_temp;
                index_e += 1;
                if(index_e == miu_ek->used){
                    index_e = 0;
                    index_delta_e += 1;
                }
            }
        }
        
    }
    return weight;
}

void regionizeMiu(ArrayFloat *miu_e, ArrayFloat *miu_delta, ArrayString *member_e, ArrayString *member_delta){
    int e_idx = 0;
    
    if(miu_e->used < 2){
        if((strcmp(member_e->array[0], "N") == 0)){
            insertArrayFloat(miu_e, 0.0);
        }
        else if((strcmp(member_e->array[0], "P") == 0)){
            float temp_miu_e;
            temp_miu_e = miu_e->array[0];
            freeArrayFloat(miu_e);
            initArrayFloat(miu_e, 1);

            insertArrayFloat(miu_e, 0.0);
            insertArrayFloat(miu_e, temp_miu_e);
        }
        else if((strcmp(member_e->array[0], "Z") == 0)){
            if(miu_delta->used < 2){
                if((strcmp(member_delta->array[0], "N") == 0)){
                    float temp_miu_e;
                    temp_miu_e = miu_e->array[0];
                    freeArrayFloat(miu_e);
                    initArrayFloat(miu_e, 1);

                    insertArrayFloat(miu_e, 0.0);
                    insertArrayFloat(miu_e, temp_miu_e);
                }
                else if((strcmp(member_delta->array[0], "P") == 0)){
                    insertArrayFloat(miu_e, 0.0);
                }
                else if((strcmp(member_delta->array[0], "Z") == 0)){
                    float temp_miu_e;
                    temp_miu_e = miu_e->array[0];
                    freeArrayFloat(miu_e);
                    initArrayFloat(miu_e, 1);

                    insertArrayFloat(miu_e, 0.0);
                    insertArrayFloat(miu_e, temp_miu_e);
                }
            }
            else if(miu_delta->used>=2){
                if((strcmp(member_delta->array[0], "N") == 0)){
                    float temp_miu_e;
                    temp_miu_e = miu_e->array[0];
                    freeArrayFloat(miu_e);
                    initArrayFloat(miu_e, 1);

                    insertArrayFloat(miu_e, 0.0);
                    insertArrayFloat(miu_e, temp_miu_e);
                }
                else if((strcmp(member_delta->array[0], "Z") == 0)){
                    float temp_miu_e;
                    temp_miu_e = miu_e->array[0];
                    freeArrayFloat(miu_e);
                    initArrayFloat(miu_e, 1);

                    insertArrayFloat(miu_e, 0.0);
                    insertArrayFloat(miu_e, temp_miu_e);
                }
            }
        }
    }
    if(miu_delta->used < 2){
        if((strcmp(member_delta->array[0], "N") == 0)){
            insertArrayFloat(miu_delta, 0.0);
        }
        else if((strcmp(member_delta->array[0], "P") == 0)){
            float temp_miu_delta;
            temp_miu_delta = miu_delta->array[0];
            freeArrayFloat(miu_delta);
            initArrayFloat(miu_delta, 1);

            insertArrayFloat(miu_delta, 0.0);
            insertArrayFloat(miu_delta, temp_miu_delta);
        }
        else if((strcmp(member_delta->array[0], "Z") == 0)){
            if(miu_e->used < 2){
                if((strcmp(member_e->array[0], "N") == 0)){
                    float temp_miu_delta;
                    temp_miu_delta = miu_delta->array[0];
                    freeArrayFloat(miu_delta);
                    initArrayFloat(miu_delta, 1);

                    insertArrayFloat(miu_delta, 0.0);
                    insertArrayFloat(miu_delta, temp_miu_delta);
                }
                else if((strcmp(member_e->array[0], "P") == 0)){
                    insertArrayFloat(miu_delta, 0.0);
                }
                else if((strcmp(member_e->array[0], "Z") == 0)){
                    float temp_miu_delta;
                    temp_miu_delta = miu_delta->array[0];
                    freeArrayFloat(miu_delta);
                    initArrayFloat(miu_delta, 1);

                    insertArrayFloat(miu_delta, 0.0);
                    insertArrayFloat(miu_delta, temp_miu_delta);
                    return;
                }
            }
            else if(miu_e->used>=2){
                if((strcmp(member_e->array[0], "N") == 0)){
                    float temp_miu_delta;
                    temp_miu_delta = miu_delta->array[0];
                    freeArrayFloat(miu_delta);
                    initArrayFloat(miu_delta, 1);

                    insertArrayFloat(miu_delta, 0.0);
                    insertArrayFloat(miu_delta, temp_miu_delta);
                }
                else if((strcmp(member_e->array[0], "Z") == 0)){
                    float temp_miu_delta;
                    temp_miu_delta = miu_delta->array[0];
                    freeArrayFloat(miu_delta);
                    initArrayFloat(miu_delta, 1);

                    insertArrayFloat(miu_delta, 0.0);
                    insertArrayFloat(miu_delta, temp_miu_delta);
                    return;
                }
                else{
                    printf("EHHHH??: %s\n", member_e->array[0]);
                }
            }
        }
    }
}
Array checkNetTrue(ArrayFloat *rule_true){
    Array region;
    initArray(&region, 1);
    if(rule_true->array[0] > 0){
        insertArray(&region, 0);
        insertArray(&region, 1);
        insertArray(&region, 3);
        insertArray(&region, 4);
    }
    else if(rule_true->array[2] > 0){
        insertArray(&region, 1);
        insertArray(&region, 2);
        insertArray(&region, 4);
        insertArray(&region, 5);
    }
    else if(rule_true->array[6] > 0){
        insertArray(&region, 3);
        insertArray(&region, 4);
        insertArray(&region, 6);
        insertArray(&region, 7);
    }
    else if(rule_true->array[8] > 0){
        insertArray(&region, 4);
        insertArray(&region, 5);
        insertArray(&region, 7);
        insertArray(&region, 8);
    }
    else if(rule_true->array[5] > 0){
        insertArray(&region, 4);
        insertArray(&region, 5);
        insertArray(&region, 7);
        insertArray(&region, 8);
    }
    else if(rule_true->array[7] > 0){
        insertArray(&region, 4);
        insertArray(&region, 5);
        insertArray(&region, 7);
        insertArray(&region, 8);
    }
    else{
        insertArray(&region, 0);
        insertArray(&region, 1);
        insertArray(&region, 3);
        insertArray(&region, 4);   
    }

    for(int i=0;i<4;i++){
        rule_true->array[region.array[i]] = 1.0;
    }

    return region;
}
ArrayFloat normalize(ArrayFloat *weight,float (*omega)[4], float R[9], ArrayFloat *net_true, Array* region) {
    ArrayFloat net;
    initArrayFloat(&net, 1);
    float net_temp;

    float total_weight = 0.0;
    // printf("used net: %d", net.used);

    for(int n=0;n<4;n++){
        total_weight = 0.0;
        for(int m=0;m<4;m++){
            total_weight += weight->array[region->array[m]] * omega[region->array[m]][n];
        }
        // printf("tot weight%d: %f\n", n, total_weight);
        net_temp = (weight->array[region->array[n]]*omega[region->array[n]][n])/total_weight;
        insertArrayFloat(&net, net_temp);
        if(net_temp != 0.0){
            net_true->array[n] = 1.0;
        }
    }

    return net;
}

float normalizeToR(ArrayFloat *net, float R[9], Array *region){
    float total = 0.0;
    for(int n=0;n<4;n++){
        total += net->array[n]*R[region->array[n]];
        
    }
    // printf("total: %f\n", total);

    return total;
}

float calculateCost(float target, float output){
    return pow((target-output), 2);
}
void printArray(Array *a) {
    printf("-------------\n");
  for (size_t i = 0; i < a->used; i++) {
    printf("%d\n", a->array[i]);
  }
  printf("-------------\n");
}
void printArrayFloat(ArrayFloat *a) {
    printf("-------------\n");
  for (size_t i = 0; i < a->used; i++) {
    printf("%f\n", a->array[i]);
  }
  printf("-------------\n");
}
void printArrayString(ArrayString *a) {
    printf("-------------\n");
  for (size_t i = 0; i < a->used; i++) {
    printf("%s\n", a->array[i]);
  }
  printf("-------------\n");
}
float x1(float prev_x1, float prev_x2, float u){
    return prev_x1 + 0.01*prev_x2 + 0.01*u;
}
    

float x2(float prev_x1, float prev_x2, float u){
    return 0.1*prev_x1 + 0.97*prev_x2;
}

float r(float k){
    if(k<500)
        return sin(M_PI*k/25);
    else if((k>=500) && (k < 1000))
        return 1;
    else if((k>=1000) && (k<1500))
        return -1;
    else if(k>=1500)
        return 0.3*sin(M_PI*k/25) + 0.4*sin(M_PI*k/32) + 0.3*sin(M_PI*k/40);
}

int main(){
    float omega[9][4];
    float R[9] = {-10.0, -10.0, 10.0, -10.0, 0.0, 10.0, -10.0, 10.0, 10.0};
    for(int m=0;m<9;m++){
        for(int n=0;n<4;n++){
            omega[m][n] = 1.0;
        }
    }
    ArrayFloat rules_true;
    initArrayFloat(&rules_true, 1);
    for(int m=0; m<9;m++){
        insertArrayFloat(&rules_true, 0.0);
    }
    ArrayFloat net_true;
    initArrayFloat(&net_true, 1);
    for(int n=0; n<4;n++){
        insertArrayFloat(&net_true, 0.0);
    }

    float output_error;

    ArrayString membership_ek;
    ArrayFloat miu_ek;

    ArrayString membership_delta_ek;
    ArrayFloat miu_delta_ek;

    float initial_y = 0.0;
    float desired_output = 1.0;

    int t[TIME];
    for(int i=0;i<TIME;i++){
        t[i] = i;
    }

    
    ArrayFloat u;
    initArrayFloat(&u, 1);

    ArrayFloat u_delta;
    initArrayFloat(&u_delta, 1);

    ArrayFloat y;
    initArrayFloat(&y, 1);
    insertArrayFloat(&y, 0.0);

    ArrayFloat y_delta;
    initArrayFloat(&y_delta, 1);
    insertArrayFloat(&y_delta, 0.0);

    ArrayFloat y_m;
    initArrayFloat(&y_m, 1);
    insertArrayFloat(&y_m, 0.0);

    ArrayFloat x_1;
    initArrayFloat(&x_1, 1);
    insertArrayFloat(&x_1, 0.0);

    ArrayFloat x_2;
    initArrayFloat(&x_2, 1);
    insertArrayFloat(&x_2, 0.0);

    float e = 0.0;
    float prev_e = 0.0;
    float delta_e = 0.0;
    float delta2_e = 0.0;
    float prev_delta_e = 0.0;

    // delta uk
    for(int epoch=0;epoch<EPOCH;epoch++){
        for(int k=0;k<TIME;k++){
            desired_output = y_m.array[k];
            if(y_m.size > 1){
                float y_m_res = 0.6*y_m.array[k] + 0.2*y_m.array[k-1] + 0.1*r(k);
                insertArrayFloat(&y_m, y_m_res);
            }
            else{
                float y_m_res = 0.3*y_m.array[k] + 0.1*r(k);
                insertArrayFloat(&y_m, y_m_res);
            }
            if(k != 0){
                prev_e = e;
                prev_delta_e = delta_e;
            }
            else{
                prev_e = desired_output - y_delta.array[k];
            }
            e = desired_output - y_delta.array[k];
            delta_e = e - prev_e;
            delta2_e = delta_e - prev_delta_e;

            // printf("desired output: %f\n", desired_output);
            // printf("error: %f\n", e);

            initArrayString(&membership_ek, 1);
            initArrayFloat(&miu_ek, 1);

            initArrayString(&membership_delta_ek, 1);
            initArrayFloat(&miu_delta_ek, 1);

            fuzzifier fuzzifier_result;
            fuzzifier_result = fuzzifierFunc(e, -1.0, 0.0, 1.0);

            // printf("u membership and miu: \n");
            // printArrayString(&fuzzifier_result.membership);
            // printArrayFloat(&fuzzifier_result.miu);
            
            for(size_t i=0; i<fuzzifier_result.membership.used;i++){
                insertArrayString(&membership_ek, fuzzifier_result.membership.array[i]);
                insertArrayFloat(&miu_ek, fuzzifier_result.miu.array[i]);
            }
            fuzzifier fuzzifier_result_delta = fuzzifierFunc(delta_e, -1.0, 0.0, 1.0);

            // printf("delta u membership and miu:\n");
            // printArrayString(&fuzzifier_result_delta.membership);
            // printArrayFloat(&fuzzifier_result_delta.miu);
            
            for(size_t i=0; i<fuzzifier_result_delta.membership.used;i++){
                insertArrayString(&membership_delta_ek, fuzzifier_result_delta.membership.array[i]);
                insertArrayFloat(&miu_delta_ek, fuzzifier_result_delta.miu.array[i]);
            }

            regionizeMiu(&miu_ek, &miu_delta_ek, &membership_ek, &membership_delta_ek);

            ArrayString fired_rules;
            ArrayFloat weight;
            ArrayFloat net;
            Array region;

            fired_rules = rule(&membership_ek, &membership_delta_ek, &rules_true);
            region = checkNetTrue(&rules_true);

            // printArrayFloat(&miu_ek);
            // printArrayFloat(&miu_delta_ek);
            // printf("rules true\n");
            // printArrayFloat(&rules_true);
            // printArrayString(&fired_rules);
            weight = defuzzification(&miu_ek, &miu_delta_ek, &rules_true);

            // printf("weight: \n");
            // printArrayFloat(&weight);

            net = normalize(&weight,omega, R, &net_true, &region);
            // printf("net:\n");
            // printArrayFloat(&net);
            // printArrayFloat(&net_true);

            float u_delta_res  = normalizeToR(&net, R, &region);
            // printf("u_delta res: %f\n", u_delta_res);
            insertArrayFloat(&u_delta, u_delta_res);

            if(u_delta.size > 1){
                u_delta.array[k] = u_delta.array[k] + u_delta.array[k-1];

                float y_delta_res = 
                        0.35*((y_delta.array[k]*y_delta.array[k-1] * (y_delta.array[k]+2.5))/
                        (1+pow(y_delta.array[k],2) + pow(y_delta.array[k-1],2))) + 0.35*u_delta.array[k];
                if(fabs(y_delta_res)>100){
                    printf("\nerror at k=%d,epoch=%d\n", k,epoch);
                    printf("desired output: %f\n", desired_output);
                    printf("error: %f\n", e);
                    printf("u membership and miu: \n");
                    printArrayFloat(&miu_ek);
                    printArrayString(&membership_ek);
                    printf("u membership and miu: \n");
                    printArrayFloat(&miu_delta_ek);
                    printArrayString(&membership_delta_ek);
                    printf("rules true\n");
                    printArrayFloat(&rules_true);
                    printArrayString(&fired_rules);
                    printf("weight: \n");
                    printArrayFloat(&weight);
                    printf("net:\n");
                    printArrayFloat(&net);
                    printArrayFloat(&net_true);
                    return 1;
                }
                insertArrayFloat(&y_delta, y_delta_res);
            }
            else{
                float y_delta_res = 
                    0.35*((y_delta.array[k]*y_delta.array[k-1] * (y_delta.array[k]+2.5))/
                    (1+pow(y_delta.array[k],2) + pow(y_delta.array[k-1],2))) + 0.35*u_delta.array[k];
                insertArrayFloat(&y_delta, y_delta_res);
            }

            
            // printf("u: %f\n", u_delta.array[k]);
            // printf("y: %f\n", y_delta.array[k+1]);
            output_error = (y_delta.array[k+1] - desired_output) * sgn((y_delta.array[k+1] - y_delta.array[k])/(u_delta.array[k] - u_delta.array[k-1]));
            

            // printf("output err: %f\n", output_error);

            float hidded_error = 0.0;
            for(int i=0;i<9;i++){
                hidded_error += output_error * R[i];
            }

            // printf("hidden err: %f\n", hidded_error);
            // printf("OK\n");

            float delta_R;
            int idx_net=0;
            for(int i=0;i<9;i++){
                if(rules_true.array[i] > 0){
                    if(idx_net>3){
                        printf("EHHHH?\n");
                        return 1;
                    }
                    if(i==0){
                        R[0] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[1] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[3] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[6] -= LEARNING_RATE*output_error*(net.array[idx_net]);

                        R[0] = fmin(R[0], -0.1);
                        R[1] = fmin(R[1], -0.1);
                        R[3] = fmin(R[3], -0.1);
                        R[6] = fmin(R[6], -0.1);
                        
                    }
                    else if(i==1){
                        R[0] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[1] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[3] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[6] -= LEARNING_RATE*output_error*(net.array[idx_net]);

                        R[0] = fmin(R[0], -0.1);
                        R[1] = fmin(R[1], -0.1);
                        R[3] = fmin(R[3], -0.1);
                        R[6] = fmin(R[6], -0.1);
                    }
                    else if(i==3){
                        R[0] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[1] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[3] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[6] -= LEARNING_RATE*output_error*(net.array[idx_net]);

                        R[0] = fmin(R[0], -0.1);
                        R[1] = fmin(R[1], -0.1);
                        R[3] = fmin(R[3], -0.1);
                        R[6] = fmin(R[6], -0.1);
                    }
                    else if(i==6){
                        R[0] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[1] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[3] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[6] -= LEARNING_RATE*output_error*(net.array[idx_net]);   

                        R[0] = fmin(R[0], -0.1);
                        R[1] = fmin(R[1], -0.1);
                        R[3] = fmin(R[3], -0.1);
                        R[6] = fmin(R[6], -0.1);
                    }

                    if(i==2){
                        R[2] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[5] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[7] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[8] -= LEARNING_RATE*output_error*(net.array[idx_net]);

                        R[2] = fmax(R[2], 0.1);
                        R[5] = fmax(R[5], 0.1);
                        R[7] = fmax(R[7], 0.1);
                        R[8] = fmax(R[8], 0.1);
                    }
                    else if(i==5){
                        R[2] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[5] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[7] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[8] -= LEARNING_RATE*output_error*(net.array[idx_net]);

                        R[2] = fmax(R[2], 0.1);
                        R[5] = fmax(R[5], 0.1);
                        R[7] = fmax(R[7], 0.1);
                        R[8] = fmax(R[8], 0.1);
                    }
                    else if(i==7){
                        R[2] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[5] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[7] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[8] -= LEARNING_RATE*output_error*(net.array[idx_net]);

                        R[2] = fmax(R[2], 0.1);
                        R[5] = fmax(R[5], 0.1);
                        R[7] = fmax(R[7], 0.1);
                        R[8] = fmax(R[8], 0.1);
                    }
                    else if(i==8){
                        R[2] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[5] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[7] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                        R[8] -= LEARNING_RATE*output_error*(net.array[idx_net]);

                        R[2] = fmax(R[2], 0.1);
                        R[5] = fmax(R[5], 0.1);
                        R[7] = fmax(R[7], 0.1);
                        R[8] = fmax(R[8], 0.1);
                    }
                    
                    if(i==4){
                        R[i] -= LEARNING_RATE*output_error*(net.array[idx_net]);
                    }
                    
                    idx_net += 1;
                }
                
            }

            float total_weight = 0.0;
            for(int i=0;i<9;i++){
                for(int n=0;n<4;n++){
                    total_weight = 0.0;
                    for(int m=0;m<9;m++){
                        total_weight += weight.array[m] * omega[m][n];
                    }
                    if(rules_true.array[i] > 0){
                        if(net_true.array[n] > 0){
                            float delta_omw = LEARNING_RATE*(output_error * R[i])*((weight.array[i]*total_weight*omega[i][n] - pow(weight.array[i], 2)*omega[i][n])/(pow(total_weight, 2))); 
                            omega[i][n] -= delta_omw;
                            // printf("OMEGA DECENT: [%d][%d]: %f\n", i, n, delta_omw);
                        }
                    }
                }
                // printf("OK");
                // printf(" \n");
            }

            // printf("COST: %f\n", calculateCost(desired_output, y_delta.array[k+1]));

            if(epoch == EPOCH-1){
                if(k<6){
                    printf("\nerror at k=%d,epoch=%d\n", k,epoch);
                    printf("desired output: %f\n", desired_output);
                    printf("error: %f\n", e);
                    printf("u membership and miu: \n");
                    printArrayFloat(&miu_ek);
                    printArrayString(&membership_ek);
                    printf("u membership and miu: \n");
                    printArrayFloat(&miu_delta_ek);
                    printArrayString(&membership_delta_ek);
                    printf("rules true\n");
                    printArrayFloat(&rules_true);
                    printArrayString(&fired_rules);
                    printf("region\n");
                    printArray(&region);
                    printf("weight: \n");
                    printArrayFloat(&weight);
                    printf("net:\n");
                    printArrayFloat(&net);
                    printArrayFloat(&net_true);
                    printf("delta u: %f\n", u_delta_res);
                    printf("u: %f\n", u_delta.array[k]);
                    printf("y: %f\n", y_delta.array[k+1]);
                }
            }

            freeArrayString(&membership_ek);
            freeArrayString(&membership_delta_ek);

            freeArrayFloat(&miu_ek);
            freeArrayFloat(&miu_delta_ek);

            freeArrayFloat(&rules_true);
            initArrayFloat(&rules_true, 1);
            for(int m=0; m<9;m++){
                insertArrayFloat(&rules_true, 0.0);
            }
            freeArrayFloat(&net_true);
            initArrayFloat(&net_true, 1);
            for(int n=0; n<4;n++){
                insertArrayFloat(&net_true, 0.0);
            }
            // printf("kos\n");
            
            
        }
        if(LEARNING_RATE > LEARNING_RATE_END){
            LEARNING_RATE *= 0.95;
        }
        
        // reset result each epoch
        if(epoch < EPOCH-1){
            freeArrayFloat(&u_delta);
            initArrayFloat(&u_delta, 1);

            freeArrayFloat(&y_delta);
            initArrayFloat(&y_delta, 1);
            insertArrayFloat(&y_delta, 0.0);

            freeArrayFloat(&y_m);
            initArrayFloat(&y_m, 1);
            insertArrayFloat(&y_m, 0.0);
        }
        // printf("SUccess\n");
    }
    

    printf("SUccess\n");
    for(int i=0;i<9;i++){
        printf("R%d: %f\n", i, R[i]);
    }
    for(int m=0;m<9;m++){
        for(int n=0;n<4;n++){
            printf("W%d%d: %f, ", m, n, omega[m][n]);
        }
        printf("\n");
    }

    double x_plot[TIME];
    double y_plot[TIME];

    double x_plot2[TIME];
    double y_plot2[TIME];

    // 0-500

    for(int i=0;i<TIME/4;i++){
        x_plot[i] = (double)t[i];
    }

    for(int i=0;i<TIME/4;i++){
        y_plot[i] = (double)y_delta.array[i];
    }

    for(int i=0;i<TIME/4;i++){
        x_plot2[i] = (double)t[i];
    }

    for(int i=0;i<TIME/4;i++){
        y_plot2[i] = (double)y_m.array[i];
    }

    RGBABitmapImageReference *imageRef = CreateRGBABitmapImageReference();
    wchar_t title[] = L"NN 3 (0 - 500)";
    DrawScatterPlot(imageRef, 600, 400, x_plot, TIME/4, y_plot, TIME/4, x_plot2, TIME/4, y_plot2, TIME/4, title);

    size_t length;
    double *pngData = ConvertToPNG(&length, imageRef->image);
    WriteToFile(pngData, length, "NN3-1.png");
    
    // 500-1000
    int idx = 0;
    int idx_plot = 500;
    while(idx<500){
        x_plot[idx] = (double)t[idx_plot];
        y_plot[idx] = (double)y_delta.array[idx_plot];
        x_plot2[idx] = (double)t[idx_plot];
        y_plot2[idx] = (double)y_m.array[idx_plot];

        idx += 1;
        idx_plot += 1;
    }

    RGBABitmapImageReference *imageRef2 = CreateRGBABitmapImageReference();
    wchar_t title2[] = L"NN 3 (500 - 1000)";
    DrawScatterPlot(imageRef2, 600, 400, x_plot, TIME/4, y_plot, TIME/4, x_plot2, TIME/4, y_plot2, TIME/4, title2);

    size_t length2;
    double *pngData2 = ConvertToPNG(&length2, imageRef2->image);
    WriteToFile(pngData2, length2, "NN3-2.png");

    // 1000-1500
    idx = 0;
    idx_plot = 1000;
    while(idx<500){
        x_plot[idx] = (double)t[idx_plot];
        y_plot[idx] = (double)y_delta.array[idx_plot];
        x_plot2[idx] = (double)t[idx_plot];
        y_plot2[idx] = (double)y_m.array[idx_plot];

        idx += 1;
        idx_plot += 1;
    }

    RGBABitmapImageReference *imageRef3 = CreateRGBABitmapImageReference();
    wchar_t title3[] = L"NN 3 (1000 - 1500)";
    DrawScatterPlot(imageRef3, 600, 400, x_plot, TIME/4, y_plot, TIME/4, x_plot2, TIME/4, y_plot2, TIME/4, title3);

    size_t length3;
    double *pngData3 = ConvertToPNG(&length3, imageRef3->image);
    WriteToFile(pngData3, length3, "NN3-3.png");

    // 1500-2000
    idx = 0;
    idx_plot = 1500;
    while(idx<500){
        x_plot[idx] = (double)t[idx_plot];
        y_plot[idx] = (double)y_delta.array[idx_plot];
        x_plot2[idx] = (double)t[idx_plot];
        y_plot2[idx] = (double)y_m.array[idx_plot];

        idx += 1;
        idx_plot += 1;
    }

    RGBABitmapImageReference *imageRef4 = CreateRGBABitmapImageReference();
    wchar_t title4[] = L"NN 3 (1500 - 2000)";
    DrawScatterPlot(imageRef4, 600, 400, x_plot, TIME/4, y_plot, TIME/4, x_plot2, TIME/4, y_plot2, TIME/4, title4);

    size_t length4;
    double *pngData4 = ConvertToPNG(&length4, imageRef4->image);
    WriteToFile(pngData4, length4, "NN3-4.png");
    
    printf("SUccess\n");
    return 0;
 }