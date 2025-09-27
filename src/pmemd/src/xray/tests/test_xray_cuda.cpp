#include <gtest/gtest.h>
#include "xray/xray_non_bulk.h"

TEST(xray_non_bulk, init_finalize) {
    pmemd_xray_non_bulk_init_gpu();
    pmemd_xray_non_bulk_finalize_gpu();
}
