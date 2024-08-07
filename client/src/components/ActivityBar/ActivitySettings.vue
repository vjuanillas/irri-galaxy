<script setup lang="ts">
import { library } from "@fortawesome/fontawesome-svg-core";
import { faSquare } from "@fortawesome/free-regular-svg-icons";
import { faCheckSquare, faThumbtack, faTrash } from "@fortawesome/free-solid-svg-icons";
import { FontAwesomeIcon } from "@fortawesome/vue-fontawesome";
import { storeToRefs } from "pinia";
import { computed, type ComputedRef } from "vue";

import { type Activity, useActivityStore } from "@/stores/activityStore";

library.add({
    faCheckSquare,
    faSquare,
    faTrash,
    faThumbtack,
});

const props = defineProps<{
    query: string;
}>();

const activityStore = useActivityStore();
const { activities } = storeToRefs(activityStore);

const filteredActivities = computed(() => {
    if (props.query?.length > 0) {
        const queryLower = props.query.toLowerCase();
        const results = activities.value.filter((a: Activity) => {
            const attributeValues = [a.title, a.description];
            for (const value of attributeValues) {
                if (value.toLowerCase().indexOf(queryLower) !== -1) {
                    return true;
                }
            }
            return false;
        });
        return results;
    } else {
        return activities.value;
    }
});

const foundActivities: ComputedRef<boolean> = computed(() => {
    return filteredActivities.value.length > 0;
});

function onClick(activity: Activity) {
    if (activity.optional) {
        activity.visible = !activity.visible;
    }
}

function onRemove(activity: Activity) {
    activityStore.remove(activity.id);
}
</script>

<template>
    <div class="activity-settings rounded no-highlight">
        <div v-if="foundActivities" class="activity-settings-content">
            <div v-for="activity in filteredActivities" :key="activity.id">
                <button class="activity-settings-item p-2 cursor-pointer" @click="onClick(activity)">
                    <div class="d-flex justify-content-between align-items-start">
                        <span class="w-100">
                            <FontAwesomeIcon
                                v-if="!activity.optional"
                                class="icon-check mr-1"
                                icon="fas fa-thumbtack"
                                fa-fw />
                            <FontAwesomeIcon
                                v-else-if="activity.visible"
                                class="icon-check mr-1"
                                icon="fas fa-check-square"
                                fa-fw />
                            <FontAwesomeIcon v-else class="mr-1" icon="far fa-square" fa-fw />
                            <span>
                                <icon class="mr-1" :icon="activity.icon" />
                                <span v-localize class="font-weight-bold">{{
                                    activity.title || "No title available"
                                }}</span>
                            </span>
                        </span>
                        <b-button
                            v-if="activity.mutable"
                            data-description="delete activity"
                            class="button-delete"
                            size="sm"
                            variant="link"
                            @click.stop="onRemove(activity)">
                            <FontAwesomeIcon icon="fa-trash" fa-fw />
                        </b-button>
                    </div>
                    <div v-localize class="text-muted">
                        {{ activity.description || "No description available" }}
                    </div>
                </button>
            </div>
        </div>
        <div v-else class="activity-settings-content">
            <b-alert v-localize class="py-1 px-2" show> No matching activities found. </b-alert>
        </div>
    </div>
</template>

<style lang="scss">
@import "theme/blue.scss";

.activity-settings {
    overflow-y: hidden;
    display: flex;
    flex-direction: column;
}

.activity-settings-content {
    overflow-y: auto;
}

.activity-settings-item {
    background: none;
    border: none;
    text-align: left;
    transition: none;
    width: 100%;

    .icon-check {
        color: darken($brand-success, 15%);
    }
    .button-delete {
        background: transparent;
    }
}
.activity-settings-item:hover {
    background: $gray-200;
}
</style>
